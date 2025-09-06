"""SVG renderer for gene review visualizations."""

import drawsvg as draw
from dataclasses import dataclass
from typing import Dict, Any, Optional, Tuple
from .layout_engine import LayoutConfig, BoundingBox, Position


@dataclass
class SVGRenderer:
    """Render gene review visualizations as SVG.
    
    This class uses drawsvg to create clean, professional SVG visualizations
    of gene annotation reviews with action indicators.
    
    Examples:
        >>> config = LayoutConfig()
        >>> renderer = SVGRenderer(config)
        >>> # Create sample data for rendering
        >>> layout_data = {
        ...     'term_positions': {},
        ...     'slim_positions': {},
        ...     'namespace_positions': {},
        ...     'action_positions': {},
        ...     'reason_positions': {},
        ...     'support_positions': {},
        ...     'canvas_dimensions': (1800, 800)
        ... }
        >>> annotations = []
        >>> gene_info = {'symbol': 'TEST', 'name': 'Test Gene'}
        >>> from unittest.mock import Mock
        >>> slim_mapper = Mock()
        >>> slim_mapper.get_term_label = Mock(return_value="Test Label")
        >>> svg = renderer.render(layout_data, annotations, gene_info, slim_mapper)
        >>> isinstance(svg, draw.Drawing)
        True
    """
    
    config: LayoutConfig
    
    def render(self, 
              layout: Dict[str, Any],
              annotations: list,
              gene_info: Dict[str, str],
              slim_mapper: Any) -> draw.Drawing:
        """Render the complete visualization.
        
        Args:
            layout: Layout information from LayoutEngine
            annotations: List of existing annotations
            gene_info: Gene metadata (symbol, taxon, etc.)
            slim_mapper: GO slim mapper for labels
            
        Returns:
            drawsvg Drawing object
        """
        # Get canvas dimensions from layout
        width, height = layout["canvas_dimensions"]
        
        # Create drawing
        d = draw.Drawing(width, height)
        
        # Add background
        d.append(draw.Rectangle(0, 0, width, height, fill=self.config.colors["background"]))
        
        # Add title
        self._render_title(d, gene_info, width)
        
        # Render namespace headers
        for namespace, box in layout["namespace_positions"].items():
            self._render_namespace_header(d, namespace, box)
        
        # Render slim categories
        for slim_id, box in layout["slim_positions"].items():
            label = slim_mapper.get_term_label(slim_id) or slim_id
            self._render_slim_category(d, label, box)
        
        # Render term boxes and actions
        annotation_map = self._create_annotation_map(annotations)
        
        for term_id, box in layout["term_positions"].items():
            # Skip new term boxes (handled with modifications)
            if term_id.endswith("_new"):
                continue
                
            # Get annotation for this term
            annotation = annotation_map.get(term_id)
            
            # Get term label
            if annotation and hasattr(annotation.term, 'label'):
                label = annotation.term.label
            else:
                label = slim_mapper.get_term_label(term_id) or term_id
            
            self._render_term_box(d, term_id, label, box, annotation)
        
        # Render actions (arrows, symbols)
        for action_id, (start, end, action_type) in layout["action_positions"].items():
            self._render_action(d, start, end, action_type, action_id, annotation_map, slim_mapper, layout)
        
        # Render reason boxes
        for term_id, box in layout.get("reason_positions", {}).items():
            annotation = annotation_map.get(term_id)
            if annotation and hasattr(annotation, 'review') and annotation.review:
                self._render_reason_box(d, box, annotation.review)
        
        # Render support boxes
        for term_id, boxes in layout.get("support_positions", {}).items():
            annotation = annotation_map.get(term_id)
            if annotation and hasattr(annotation, 'review') and annotation.review:
                if hasattr(annotation.review, 'supported_by') and annotation.review.supported_by:
                    for i, box in enumerate(boxes):
                        if i < len(annotation.review.supported_by):
                            self._render_support_box(d, box, annotation.review.supported_by[i])
        
        # Render legend
        legend_box = BoundingBox(
            x=(width - 600) / 2,
            y=height - self.config.margin_bottom - 60,
            width=600,
            height=60
        )
        self._render_legend(d, legend_box)
        
        return d
    
    def _create_annotation_map(self, annotations: list) -> Dict[str, Any]:
        """Create a map of term IDs to annotations."""
        annotation_map = {}
        for annotation in annotations:
            if hasattr(annotation, 'term') and hasattr(annotation.term, 'id'):
                annotation_map[annotation.term.id] = annotation
        return annotation_map
    
    def _render_title(self, d: draw.Drawing, gene_info: Dict[str, str], width: float):
        """Render the title section."""
        title_text = f"Gene Review: {gene_info.get('gene_symbol', 'Unknown')}"
        if 'taxon' in gene_info and 'label' in gene_info['taxon']:
            title_text += f" ({gene_info['taxon']['label']})"
        
        d.append(draw.Text(
            title_text,
            x=width / 2,
            y=30,
            font_size=self.config.title_font_size,
            text_anchor="middle",
            font_family="Arial, sans-serif",
            font_weight="bold",
            fill=self.config.colors["text_color"]
        ))
    
    def _render_namespace_header(self, d: draw.Drawing, namespace: str, box: BoundingBox):
        """Render a namespace header."""
        # Background rectangle
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill=self.config.colors["namespace_bg"],
            stroke=self.config.colors["border_color"],
            stroke_width=1,
            rx=5, ry=5
        ))
        
        # Label
        label = namespace.replace("_", " ").title()
        d.append(draw.Text(
            label,
            x=box.x + 10,
            y=box.y + box.height / 2 + 5,
            font_size=self.config.namespace_font_size,
            font_family="Arial, sans-serif",
            font_weight="bold",
            fill=self.config.colors["text_color"]
        ))
    
    def _render_slim_category(self, d: draw.Drawing, label: str, box: BoundingBox):
        """Render a GO slim category box."""
        # Background rectangle
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill=self.config.colors["slim_bg"],
            stroke=self.config.colors["border_color"],
            stroke_width=1,
            rx=3, ry=3
        ))
        
        # Label - increased limit for longer category names
        d.append(draw.Text(
            self._truncate_text(label, 40),  # Further increased to use available space
            x=box.x + 8,
            y=box.y + box.height / 2 + 4,
            font_size=self.config.slim_font_size,
            font_family="Arial, sans-serif",
            fill=self.config.colors["text_color"]
        ))
    
    def _render_term_box(self, d: draw.Drawing, term_id: str, label: str, 
                        box: BoundingBox, annotation: Optional[Any] = None):
        """Render a GO term box."""
        # Determine box color based on annotation status
        fill_color = self.config.colors["term_bg"]
        border_color = self.config.colors["border_color"]
        border_width = 1
        
        if annotation and hasattr(annotation, 'review') and annotation.review:
            action = annotation.review.action if hasattr(annotation.review, 'action') else None
            if action:
                # Adjust border color based on action
                action_colors = {
                    "ACCEPT": self.config.colors["accept_color"],
                    "REMOVE": self.config.colors["remove_color"],
                    "MODIFY": self.config.colors["modify_color"],
                    "MARK_AS_OVER_ANNOTATED": self.config.colors["over_annotated_color"],
                    "UNDECIDED": self.config.colors["undecided_color"],
                    "KEEP_AS_NON_CORE": self.config.colors["non_core_color"]
                }
                border_color = action_colors.get(str(action), border_color)
                border_width = 2
        
        # Background rectangle
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill=fill_color,
            stroke=border_color,
            stroke_width=border_width,
            rx=3, ry=3
        ))
        
        # Term ID
        d.append(draw.Text(
            term_id,
            x=box.x + 5,
            y=box.y + 14,
            font_size=self.config.label_font_size,
            font_family="monospace",
            fill="#666666"
        ))
        
        # Term label - increased character limit
        d.append(draw.Text(
            self._truncate_text(label, 38),  # Much more space available in 200px width
            x=box.x + 5,
            y=box.y + 30,
            font_size=self.config.term_font_size,
            font_family="Arial, sans-serif",
            fill=self.config.colors["text_color"]
        ))
    
    def _render_action(self, d: draw.Drawing, start: Position, end: Position, 
                      action_type: str, action_id: str, annotation_map: Dict, slim_mapper: Any,
                      layout: Dict[str, Any]):
        """Render an action indicator (arrow, symbol, etc.)."""
        if action_type == "accept":
            # Green checkmark
            self._render_checkmark(d, end.x - 20, end.y, self.config.colors["accept_color"])
            
        elif action_type == "remove":
            # Red X
            self._render_x_mark(d, end.x - 20, end.y, self.config.colors["remove_color"])
            
        elif action_type == "modify":
            # Arrow to new term
            self._render_arrow(d, start, end, self.config.colors["modify_color"])
            
            # Extract the original term ID from the action_id (format: "GO:xxxxx_modify")
            if action_id.endswith("_modify"):
                base_term_id = action_id[:-7]  # Remove "_modify" to get original term ID
                
                # Get the annotation to find the proposed replacement terms
                annotation = annotation_map.get(base_term_id)
                if annotation and hasattr(annotation, 'review') and annotation.review:
                    review = annotation.review
                    if hasattr(review, 'proposed_replacement_terms') and review.proposed_replacement_terms:
                        # Get the first replacement term
                        new_term = review.proposed_replacement_terms[0]
                        new_term_id = new_term.id if hasattr(new_term, 'id') else None
                        
                        if new_term_id:
                            # Look for the new term box using the unique key
                            unique_key = f"{base_term_id}_to_{new_term_id}_new"
                            if unique_key in layout["term_positions"]:
                                box = layout["term_positions"][unique_key]
                                label = new_term.label if hasattr(new_term, 'label') else (slim_mapper.get_term_label(new_term_id) or new_term_id)
                                
                                # Render new term box with different background
                                d.append(draw.Rectangle(
                                    box.x, box.y, box.width, box.height,
                                    fill=self.config.colors["new_term_bg"],
                                    stroke=self.config.colors["modify_color"],
                                    stroke_width=2,
                                    rx=3, ry=3
                                ))
                                
                                # Term ID
                                d.append(draw.Text(
                                    new_term_id,
                                    x=box.x + 5,
                                    y=box.y + 14,
                                    font_size=self.config.label_font_size,
                                    font_family="monospace",
                                    fill="#666666"
                                ))
                                
                                # Term label - use same increased limit as regular term boxes
                                d.append(draw.Text(
                                    self._truncate_text(label, 38),  # Consistent with regular term boxes
                                    x=box.x + 5,
                                    y=box.y + 30,
                                    font_size=self.config.term_font_size,
                                    font_family="Arial, sans-serif",
                                    fill=self.config.colors["text_color"]
                                ))
            
        elif action_type == "over_annotated":
            # Warning symbol
            self._render_warning(d, end.x - 20, end.y, self.config.colors["over_annotated_color"])
            
        elif action_type == "non_core":
            # Hollow circle
            self._render_hollow_circle(d, end.x - 20, end.y, self.config.colors["non_core_color"])
            
        elif action_type == "undecided":
            # Question mark
            self._render_question_mark(d, end.x - 20, end.y, self.config.colors["undecided_color"])
    
    def _render_checkmark(self, d: draw.Drawing, x: float, y: float, color: str):
        """Render a checkmark symbol."""
        path = draw.Path(stroke=color, stroke_width=3, fill="none")
        path.M(x - 8, y)
        path.L(x - 3, y + 5)
        path.L(x + 8, y - 8)
        d.append(path)
    
    def _render_x_mark(self, d: draw.Drawing, x: float, y: float, color: str):
        """Render an X mark."""
        d.append(draw.Line(x - 8, y - 8, x + 8, y + 8, stroke=color, stroke_width=3))
        d.append(draw.Line(x - 8, y + 8, x + 8, y - 8, stroke=color, stroke_width=3))
    
    def _render_arrow(self, d: draw.Drawing, start: Position, end: Position, color: str):
        """Render an arrow."""
        # Arrow line
        d.append(draw.Line(start.x, start.y, end.x, end.y, stroke=color, stroke_width=2))
        
        # Arrowhead
        path = draw.Path(fill=color)
        path.M(end.x, end.y)
        path.L(end.x - 8, end.y - 4)
        path.L(end.x - 8, end.y + 4)
        path.Z()
        d.append(path)
    
    def _render_warning(self, d: draw.Drawing, x: float, y: float, color: str):
        """Render a warning triangle."""
        path = draw.Path(fill=color)
        path.M(x, y - 10)
        path.L(x - 10, y + 8)
        path.L(x + 10, y + 8)
        path.Z()
        d.append(path)
        
        # Exclamation mark
        d.append(draw.Text(
            "!",
            x=x,
            y=y + 4,
            font_size=14,
            text_anchor="middle",
            font_family="Arial, sans-serif",
            font_weight="bold",
            fill="white"
        ))
    
    def _render_hollow_circle(self, d: draw.Drawing, x: float, y: float, color: str):
        """Render a hollow circle."""
        d.append(draw.Circle(x, y, 8, fill="none", stroke=color, stroke_width=2))
    
    def _render_question_mark(self, d: draw.Drawing, x: float, y: float, color: str):
        """Render a question mark."""
        d.append(draw.Text(
            "?",
            x=x,
            y=y + 5,
            font_size=20,
            text_anchor="middle",
            font_family="Arial, sans-serif",
            font_weight="bold",
            fill=color
        ))
    
    def _render_legend(self, d: draw.Drawing, box: BoundingBox):
        """Render the legend."""
        # Background
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill="white",
            stroke=self.config.colors["border_color"],
            stroke_width=1,
            rx=5, ry=5
        ))
        
        # Legend items
        items = [
            ("✓", "Accept", self.config.colors["accept_color"]),
            ("→", "Modify", self.config.colors["modify_color"]),
            ("✗", "Remove", self.config.colors["remove_color"]),
            ("⚠", "Over-annotated", self.config.colors["over_annotated_color"]),
            ("○", "Non-core", self.config.colors["non_core_color"]),
            ("?", "Undecided", self.config.colors["undecided_color"])
        ]
        
        x_offset = box.x + 20
        y_center = box.y + box.height / 2
        
        for i, (symbol, label, color) in enumerate(items):
            x = x_offset + i * 95
            
            # Symbol
            d.append(draw.Text(
                symbol,
                x=x,
                y=y_center + 5,
                font_size=16,
                font_family="Arial, sans-serif",
                fill=color
            ))
            
            # Label
            d.append(draw.Text(
                label,
                x=x + 20,
                y=y_center + 5,
                font_size=12,
                font_family="Arial, sans-serif",
                fill=self.config.colors["text_color"]
            ))
    
    def _truncate_text(self, text: str, max_length: int) -> str:
        """Truncate text to fit in boxes."""
        if len(text) <= max_length:
            return text
        return text[:max_length - 3] + "..."
    
    def _render_reason_box(self, d: draw.Drawing, box: BoundingBox, review: Any):
        """Render a reason text box."""
        # Background with light yellow color
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill="#fff9c4",  # Light yellow
            stroke="#f9a825",  # Amber border
            stroke_width=1,
            rx=3, ry=3
        ))
        
        # Header
        d.append(draw.Text(
            "Reason:",
            x=box.x + 5,
            y=box.y + 10,
            font_size=self.config.reason_font_size,
            font_family="Arial, sans-serif",
            font_weight="bold",
            fill="#666666"
        ))
        
        # Reason text (wrapped)
        if hasattr(review, 'reason') and review.reason:
            reason_text = str(review.reason)
            lines = self._wrap_text(reason_text, 58)  # 280px width can fit more text
            
            for i, line in enumerate(lines[:6]):  # Increased to 6 lines to use full height
                d.append(draw.Text(
                    line,
                    x=box.x + 5,
                    y=box.y + 20 + i * 9,  # Slightly tighter line spacing
                    font_size=self.config.reason_font_size,
                    font_family="Arial, sans-serif",
                    fill=self.config.colors["text_color"]
                ))
    
    def _render_support_box(self, d: draw.Drawing, box: BoundingBox, support: Any):
        """Render a supporting evidence box."""
        # Background with light blue color
        d.append(draw.Rectangle(
            box.x, box.y, box.width, box.height,
            fill="#e1f5fe",  # Light blue
            stroke="#0288d1",  # Blue border
            stroke_width=1,
            rx=3, ry=3
        ))
        
        # Extract reference info
        ref_id = ""
        ref_text = ""
        
        if hasattr(support, 'reference_id'):
            ref_id = str(support.reference_id)
        if hasattr(support, 'supporting_text'):
            ref_text = str(support.supporting_text)
        
        # Reference header
        if ref_id:
            # Make PMID shorter for display
            display_ref = ref_id.replace("PMID:", "").replace("GO_REF:", "GO:")
            if len(display_ref) > 20:
                display_ref = display_ref[:20] + "..."
                
            d.append(draw.Text(
                f"[{display_ref}]",
                x=box.x + 5,
                y=box.y + 10,
                font_size=self.config.support_font_size,
                font_family="monospace",
                font_weight="bold",
                fill="#0288d1"
            ))
        
        # Supporting text (wrapped)
        if ref_text:
            lines = self._wrap_text(ref_text, 68)  # 320px width can fit much more text
            
            start_y = box.y + 18 if ref_id else box.y + 10
            for i, line in enumerate(lines[:8]):  # Increased to 8 lines for 85px height
                d.append(draw.Text(
                    line,
                    x=box.x + 5,
                    y=start_y + i * 9,  # Tighter line spacing
                    font_size=self.config.support_font_size,
                    font_family="Arial, sans-serif",
                    fill=self.config.colors["text_color"]
                ))
    
    def _wrap_text(self, text: str, max_chars: int) -> list[str]:
        """Wrap text into lines of approximately max_chars characters."""
        words = text.split()
        lines = []
        current_line = []
        current_length = 0
        
        for word in words:
            word_length = len(word)
            if current_length + word_length + len(current_line) > max_chars:
                if current_line:
                    lines.append(" ".join(current_line))
                    current_line = [word]
                    current_length = word_length
                else:
                    # Word is too long, truncate it
                    lines.append(word[:max_chars - 3] + "...")
                    current_line = []
                    current_length = 0
            else:
                current_line.append(word)
                current_length += word_length
        
        if current_line:
            lines.append(" ".join(current_line))
        
        return lines
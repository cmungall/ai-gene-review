"""Layout engine for calculating positions of visual elements."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from ai_gene_review.datamodel.gene_review_model import ActionEnum


@dataclass
class Position:
    """2D position in SVG coordinate space."""
    x: float
    y: float


@dataclass
class BoundingBox:
    """Bounding box for a visual element."""
    x: float
    y: float
    width: float
    height: float
    
    @property
    def center(self) -> Position:
        """Get the center position of the box."""
        return Position(self.x + self.width / 2, self.y + self.height / 2)
    
    @property
    def right(self) -> float:
        """Get the right edge x-coordinate."""
        return self.x + self.width
    
    @property
    def bottom(self) -> float:
        """Get the bottom edge y-coordinate."""
        return self.y + self.height


@dataclass
class LayoutConfig:
    """Configuration for layout calculations."""
    # Canvas dimensions
    canvas_width: float = 1800  # Increased to accommodate reason and support boxes
    canvas_height: float = 800  # Will expand as needed
    
    # Margins and padding
    margin_top: float = 60
    margin_bottom: float = 80
    margin_left: float = 40
    margin_right: float = 40
    
    # Box dimensions
    term_box_width: float = 200  # Increased for longer labels
    term_box_height: float = 45  # Slightly taller for better text fit
    slim_box_width: float = 220
    slim_box_height: float = 35
    reason_box_width: float = 280  # Wider for more text
    reason_box_height: float = 70  # Taller for more lines
    support_box_width: float = 320  # Wider for supporting evidence
    support_box_height: float = 85  # Taller for more text
    
    # Spacing
    horizontal_spacing: float = 110  # Between original and modified terms
    reason_spacing: float = 35  # Between action and reason box
    support_spacing: float = 25  # Between reason and support box
    vertical_spacing: float = 100  # Much larger vertical spacing to prevent overlap
    category_spacing: float = 40  # Between different slim categories
    namespace_spacing: float = 70  # Between namespaces
    
    # Font sizes
    title_font_size: float = 20
    namespace_font_size: float = 16
    slim_font_size: float = 13
    term_font_size: float = 11
    label_font_size: float = 10
    reason_font_size: float = 9  # Smaller text for reason boxes
    support_font_size: float = 8  # Even smaller for support boxes
    
    # Colors
    colors: Dict[str, str] = field(default_factory=lambda: {
        "background": "#f5f5f5",
        "namespace_bg": "#e3f2fd",
        "slim_bg": "#e8f5e9",
        "term_bg": "#ffffff",
        "accept_color": "#4caf50",
        "modify_color": "#ff9800",
        "remove_color": "#f44336",
        "over_annotated_color": "#ffc107",
        "undecided_color": "#9e9e9e",
        "non_core_color": "#795548",
        "new_term_bg": "#c8e6c9",
        "text_color": "#212121",
        "border_color": "#bdbdbd"
    })


@dataclass
class LayoutEngine:
    """Calculate positions for visualization elements.
    
    This engine handles the layout of GO terms organized by namespace
    and slim categories, with arrows showing annotation actions.
    
    Examples:
        >>> config = LayoutConfig()
        >>> engine = LayoutEngine(config)
        >>> # Calculate layout with sample data
        >>> organized_terms = {"biological_process": {}}
        >>> annotations = []
        >>> from unittest.mock import Mock
        >>> slim_mapper = Mock()
        >>> slim_mapper.get_term_label = Mock(return_value="Test Label")
        >>> layout = engine.calculate_layout(organized_terms, annotations, slim_mapper)
        >>> isinstance(layout, dict)
        True
        >>> 'term_positions' in layout
        True
    """
    
    config: LayoutConfig = field(default_factory=LayoutConfig)
    
    # Calculated positions
    _term_positions: Dict[str, BoundingBox] = field(default_factory=dict, init=False)
    _slim_positions: Dict[str, BoundingBox] = field(default_factory=dict, init=False)
    _namespace_positions: Dict[str, BoundingBox] = field(default_factory=dict, init=False)
    _action_positions: Dict[str, Tuple[Position, Position, str]] = field(default_factory=dict, init=False)
    _reason_positions: Dict[str, BoundingBox] = field(default_factory=dict, init=False)
    _support_positions: Dict[str, List[BoundingBox]] = field(default_factory=dict, init=False)
    _canvas_height: float = field(default=800, init=False)
    _canvas_width: float = field(default=1800, init=False)
    
    def calculate_layout(self, 
                        organized_terms: Dict[str, Dict[str, List[str]]],
                        annotations: List[Any],
                        slim_mapper: Any) -> Dict[str, Any]:
        """Calculate complete layout for visualization.
        
        Args:
            organized_terms: Terms organized by namespace and slim category
            annotations: List of existing annotations with review actions
            slim_mapper: GO slim mapper for getting labels
            
        Returns:
            Dict with layout information including positions and dimensions
        """
        self._reset_positions()
        
        current_y = self.config.margin_top + 40  # Space for title
        
        # Process each namespace
        for namespace in ["biological_process", "molecular_function", "cellular_component"]:
            if namespace not in organized_terms or not organized_terms[namespace]:
                continue
                
            # Namespace header
            namespace_box = self._layout_namespace_header(namespace, current_y)
            self._namespace_positions[namespace] = namespace_box
            current_y = namespace_box.bottom + 20
            
            # Process slim categories within namespace
            for slim_id, term_ids in organized_terms[namespace].items():
                if not term_ids:
                    continue
                    
                # Slim category box
                slim_label = slim_mapper.get_term_label(slim_id) or slim_id
                slim_box = self._layout_slim_category(slim_id, slim_label, current_y)
                self._slim_positions[slim_id] = slim_box
                current_y = slim_box.bottom + 10
                
                # Layout individual terms
                for term_id in term_ids:
                    term_box = self._layout_term(term_id, current_y)
                    self._term_positions[term_id] = term_box
                    current_y = term_box.bottom + self.config.vertical_spacing
                
                current_y += self.config.category_spacing
            
            current_y += self.config.namespace_spacing
        
        # Calculate action arrows/symbols
        self._layout_actions(annotations)
        
        # Calculate maximum width based on rightmost elements
        max_width = self.config.canvas_width  # Start with default
        
        # Check all positioned elements for maximum x + width
        for positions_dict in [self._term_positions, self._reason_positions]:
            for box in positions_dict.values():
                if isinstance(box, BoundingBox):
                    max_width = max(max_width, box.right + self.config.margin_right)
        
        # Check support positions (which is a dict of lists)
        for support_list in self._support_positions.values():
            for box in support_list:
                if isinstance(box, BoundingBox):
                    max_width = max(max_width, box.right + self.config.margin_right)
        
        # Update canvas dimensions based on content
        self._canvas_height = current_y + self.config.margin_bottom + 100  # Space for legend
        self._canvas_width = max_width
        
        return {
            "term_positions": self._term_positions,
            "slim_positions": self._slim_positions,
            "namespace_positions": self._namespace_positions,
            "action_positions": self._action_positions,
            "reason_positions": self._reason_positions,
            "support_positions": self._support_positions,
            "canvas_dimensions": (self._canvas_width, self._canvas_height)
        }
    
    def _reset_positions(self):
        """Reset all position dictionaries."""
        self._term_positions = {}
        self._slim_positions = {}
        self._namespace_positions = {}
        self._action_positions = {}
        self._reason_positions = {}
        self._support_positions = {}
    
    def _layout_namespace_header(self, namespace: str, y: float) -> BoundingBox:
        """Layout a namespace header."""
        label = namespace.replace("_", " ").title()
        width = self.config.canvas_width - self.config.margin_left - self.config.margin_right
        height = 35
        
        return BoundingBox(
            x=self.config.margin_left,
            y=y,
            width=width,
            height=height
        )
    
    def _layout_slim_category(self, slim_id: str, label: str, y: float) -> BoundingBox:
        """Layout a GO slim category box."""
        return BoundingBox(
            x=self.config.margin_left + 20,  # Indent from namespace
            y=y,
            width=self.config.slim_box_width,
            height=self.config.slim_box_height
        )
    
    def _layout_term(self, term_id: str, y: float) -> BoundingBox:
        """Layout an individual GO term box."""
        return BoundingBox(
            x=self.config.margin_left + 40,  # Indent from slim category
            y=y,
            width=self.config.term_box_width,
            height=self.config.term_box_height
        )
    
    def _layout_actions(self, annotations: List[Any]):
        """Calculate positions for action arrows, reason boxes, and support boxes.
        
        Args:
            annotations: List of existing annotations with review actions
        """
        for annotation in annotations:
            if not hasattr(annotation, 'term') or not hasattr(annotation, 'review'):
                continue
                
            term_id = annotation.term.id if hasattr(annotation.term, 'id') else None
            if not term_id or term_id not in self._term_positions:
                continue
                
            review = annotation.review
            if not review or not hasattr(review, 'action'):
                continue
                
            action = review.action
            source_box = self._term_positions[term_id]
            
            # Base X position for elements to the right
            action_x = source_box.right + 80
            
            # Calculate action position based on action type
            if action == ActionEnum.ACCEPT:
                # Checkmark to the right of the term
                start = Position(source_box.right + 20, source_box.center.y)
                end = Position(source_box.right + 80, source_box.center.y)
                self._action_positions[f"{term_id}_action"] = (start, end, "accept")
                
            elif action == ActionEnum.REMOVE:
                # X mark to the right of the term
                start = Position(source_box.right + 20, source_box.center.y)
                end = Position(source_box.right + 80, source_box.center.y)
                self._action_positions[f"{term_id}_action"] = (start, end, "remove")
                
            elif action == ActionEnum.MODIFY:
                # Arrow to new term box
                if hasattr(review, 'proposed_replacement_terms') and review.proposed_replacement_terms:
                    # Create box for new term
                    new_term = review.proposed_replacement_terms[0]
                    new_term_id = new_term.id if hasattr(new_term, 'id') else None
                    if new_term_id:
                        new_box = BoundingBox(
                            x=source_box.right + self.config.horizontal_spacing,
                            y=source_box.y,
                            width=self.config.term_box_width,
                            height=self.config.term_box_height
                        )
                        # Use unique key combining original and new term IDs
                        unique_key = f"{term_id}_to_{new_term_id}_new"
                        self._term_positions[unique_key] = new_box
                        action_x = new_box.right + 20  # Position reason after new term
                        
                        # Arrow from old to new
                        start = Position(source_box.right, source_box.center.y)
                        end = Position(new_box.x, new_box.center.y)
                        self._action_positions[f"{term_id}_modify"] = (start, end, "modify")
                
            elif action == ActionEnum.MARK_AS_OVER_ANNOTATED:
                # Warning symbol to the right
                start = Position(source_box.right + 20, source_box.center.y)
                end = Position(source_box.right + 80, source_box.center.y)
                self._action_positions[f"{term_id}_action"] = (start, end, "over_annotated")
                
            elif action == ActionEnum.KEEP_AS_NON_CORE:
                # Hollow circle to the right
                start = Position(source_box.right + 20, source_box.center.y)
                end = Position(source_box.right + 80, source_box.center.y)
                self._action_positions[f"{term_id}_action"] = (start, end, "non_core")
                
            elif action == ActionEnum.UNDECIDED:
                # Question mark to the right
                start = Position(source_box.right + 20, source_box.center.y)
                end = Position(source_box.right + 80, source_box.center.y)
                self._action_positions[f"{term_id}_action"] = (start, end, "undecided")
            
            # Add reason box if reason exists
            if hasattr(review, 'reason') and review.reason:
                # Center reason box vertically with the term box
                reason_y = source_box.y + (source_box.height - self.config.reason_box_height) / 2
                reason_box = BoundingBox(
                    x=action_x + self.config.reason_spacing,
                    y=reason_y,
                    width=self.config.reason_box_width,
                    height=self.config.reason_box_height
                )
                self._reason_positions[term_id] = reason_box
                action_x = reason_box.right  # Update x position for next element
            
            # Add support boxes if supporting evidence exists
            if hasattr(review, 'supported_by') and review.supported_by:
                support_boxes = []
                # Flow support boxes horizontally instead of stacking vertically
                current_x = action_x + self.config.support_spacing
                
                for i, support in enumerate(review.supported_by):  # No limit, let them flow horizontally
                    # For multiple support boxes, place them side by side with small gap
                    support_box = BoundingBox(
                        x=current_x,
                        y=source_box.y,  # Keep aligned with the term box vertically
                        width=self.config.support_box_width,
                        height=self.config.support_box_height
                    )
                    support_boxes.append(support_box)
                    # Move x position for next support box
                    current_x = support_box.right + 10  # Small gap between support boxes
                
                if support_boxes:
                    self._support_positions[term_id] = support_boxes
    
    def get_legend_position(self) -> BoundingBox:
        """Get position for the legend."""
        legend_width = 600
        legend_height = 60
        
        return BoundingBox(
            x=(self._canvas_width - legend_width) / 2,
            y=self._canvas_height - self.config.margin_bottom - legend_height,
            width=legend_width,
            height=legend_height
        )
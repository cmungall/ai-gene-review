"""Main visualizer for gene review data."""

import yaml
from pathlib import Path
from typing import Union, Optional, Dict, Any, List
from dataclasses import dataclass
import drawsvg as draw

from ai_gene_review.datamodel.gene_review_model import GeneReview, ExistingAnnotation
from .go_slim_mapper import GOSlimMapper
from .layout_engine import LayoutEngine, LayoutConfig
from .svg_renderer import SVGRenderer


@dataclass
class ReviewVisualizer:
    """Visualize gene review annotations and their actions.
    
    This class coordinates the visualization pipeline, loading gene review data,
    organizing terms by GO slim categories, calculating layout, and rendering SVG.
    
    Examples:
        >>> visualizer = ReviewVisualizer()
        >>> # visualize_file returns a Drawing object
        >>> import os
        >>> test_file = "genes/human/CFAP300/CFAP300-ai-review.yaml"
        >>> if os.path.exists(test_file):
        ...     drawing = visualizer.visualize_file(test_file)
        ...     isinstance(drawing, draw.Drawing)
        ... else:
        ...     True  # Skip if file doesn't exist
        True
    """
    
    # Configuration
    layout_config: LayoutConfig = None
    slim_subset: str = "goslim_generic"
    
    # Components
    _slim_mapper: GOSlimMapper = None
    _layout_engine: LayoutEngine = None
    _svg_renderer: SVGRenderer = None
    _last_drawing: Optional[draw.Drawing] = None
    
    def __post_init__(self):
        """Initialize components."""
        if self.layout_config is None:
            self.layout_config = LayoutConfig()
        
        self._slim_mapper = GOSlimMapper(slim_subset=self.slim_subset)
        self._layout_engine = LayoutEngine(config=self.layout_config)
        self._svg_renderer = SVGRenderer(config=self.layout_config)
    
    def visualize_file(self, file_path: Union[str, Path]) -> draw.Drawing:
        """Visualize a gene review from a YAML file.
        
        Args:
            file_path: Path to the gene review YAML file
            
        Returns:
            drawsvg Drawing object
        """
        file_path = Path(file_path)
        
        # Load the gene review
        with open(file_path, "r") as f:
            data = yaml.safe_load(f)
        
        gene_review = GeneReview.model_validate(data)
        
        return self.visualize(gene_review)
    
    def visualize(self, gene_review: GeneReview) -> draw.Drawing:
        """Visualize a gene review object.
        
        Args:
            gene_review: GeneReview object to visualize
            
        Returns:
            drawsvg Drawing object
        """
        # Extract term IDs from existing annotations
        term_ids = self._extract_term_ids(gene_review.existing_annotations)
        
        # Organize terms by GO slim categories
        organized_terms = self._slim_mapper.organize_by_slim(term_ids)
        
        # Calculate layout
        layout = self._layout_engine.calculate_layout(
            organized_terms,
            gene_review.existing_annotations,
            self._slim_mapper
        )
        
        # Prepare gene info for rendering
        gene_info = {
            "gene_symbol": gene_review.gene_symbol,
            "taxon": {
                "id": gene_review.taxon.id if gene_review.taxon else None,
                "label": gene_review.taxon.label if gene_review.taxon else None
            } if gene_review.taxon else {}
        }
        
        # Render SVG
        drawing = self._svg_renderer.render(
            layout,
            gene_review.existing_annotations,
            gene_info,
            self._slim_mapper
        )
        
        self._last_drawing = drawing
        return drawing
    
    def visualize_comparison(self, 
                           gene_reviews: List[GeneReview],
                           mode: str = "side_by_side") -> draw.Drawing:
        """Visualize multiple gene reviews for comparison.
        
        Args:
            gene_reviews: List of GeneReview objects to compare
            mode: Comparison mode ("side_by_side" or "overlay")
            
        Returns:
            drawsvg Drawing object
        """
        if mode == "side_by_side":
            return self._visualize_side_by_side(gene_reviews)
        elif mode == "overlay":
            return self._visualize_overlay(gene_reviews)
        else:
            raise ValueError(f"Unknown comparison mode: {mode}")
    
    def _visualize_side_by_side(self, gene_reviews: List[GeneReview]) -> draw.Drawing:
        """Create side-by-side comparison of multiple gene reviews."""
        # This would create multiple panels side by side
        # For now, just visualize the first one
        if gene_reviews:
            return self.visualize(gene_reviews[0])
        else:
            raise ValueError("No gene reviews provided")
    
    def _visualize_overlay(self, gene_reviews: List[GeneReview]) -> draw.Drawing:
        """Create overlay comparison showing differences."""
        # This would overlay multiple reviews to show differences
        # For now, just visualize the first one
        if gene_reviews:
            return self.visualize(gene_reviews[0])
        else:
            raise ValueError("No gene reviews provided")
    
    def _extract_term_ids(self, annotations: List[ExistingAnnotation]) -> List[str]:
        """Extract GO term IDs from annotations.
        
        Args:
            annotations: List of existing annotations
            
        Returns:
            List of GO term IDs
        """
        term_ids = []
        
        if not annotations:
            return term_ids
        
        for annotation in annotations:
            if hasattr(annotation, 'term') and annotation.term:
                if hasattr(annotation.term, 'id') and annotation.term.id:
                    # Only include GO terms
                    if annotation.term.id.startswith("GO:"):
                        term_ids.append(annotation.term.id)
        
        return term_ids
    
    def save(self, output_path: Union[str, Path], format: str = "svg"):
        """Save the last visualization to a file.
        
        Args:
            output_path: Path to save the visualization
            format: Output format ("svg" or "png" if Cairo is available)
        """
        if self._last_drawing is None:
            raise ValueError("No visualization to save. Call visualize() first.")
        
        output_path = Path(output_path)
        
        if format == "svg":
            self._last_drawing.save_svg(str(output_path))
        elif format == "png":
            # Requires Cairo to be installed
            try:
                self._last_drawing.save_png(str(output_path))
            except ImportError:
                raise ImportError("Cairo is required for PNG export. Install with: pip install pycairo")
        else:
            raise ValueError(f"Unknown format: {format}. Use 'svg' or 'png'")
    
    def get_summary_stats(self, gene_review: GeneReview) -> Dict[str, Any]:
        """Get summary statistics for a gene review.
        
        Args:
            gene_review: GeneReview object
            
        Returns:
            Dictionary with summary statistics
        """
        stats = {
            "total_annotations": len(gene_review.existing_annotations) if gene_review.existing_annotations else 0,
            "actions": {
                "ACCEPT": 0,
                "MODIFY": 0,
                "REMOVE": 0,
                "MARK_AS_OVER_ANNOTATED": 0,
                "KEEP_AS_NON_CORE": 0,
                "UNDECIDED": 0
            }
        }
        
        if gene_review.existing_annotations:
            for annotation in gene_review.existing_annotations:
                if hasattr(annotation, 'review') and annotation.review:
                    if hasattr(annotation.review, 'action') and annotation.review.action:
                        action_str = str(annotation.review.action)
                        if action_str in stats["actions"]:
                            stats["actions"][action_str] += 1
        
        # Calculate percentages
        total = stats["total_annotations"]
        if total > 0:
            stats["action_percentages"] = {
                action: (count / total) * 100
                for action, count in stats["actions"].items()
            }
        else:
            stats["action_percentages"] = {action: 0 for action in stats["actions"]}
        
        return stats
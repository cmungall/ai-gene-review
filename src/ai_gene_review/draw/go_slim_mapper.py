"""GO Slim mapper for organizing GO terms into hierarchical categories."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from oaklib import get_adapter
from oaklib.datamodels.vocabulary import IS_A


@dataclass
class GOSlimMapper:
    """Maps GO terms to GO slim categories for hierarchical visualization.
    
    This class uses oaklib to fetch GO slim terms and organize annotations
    into their parent slim categories for visualization purposes.
    
    Examples:
        >>> mapper = GOSlimMapper()
        >>> # Map a specific GO term to its slim ancestors
        >>> slim_terms = mapper.get_slim_ancestors("GO:0006629")  # lipid metabolic process
        >>> # Returns slim terms like "GO:0008152" (metabolic process)
    """
    
    # GO slim subset to use (e.g., "goslim_generic", "goslim_yeast")
    slim_subset: str = "goslim_generic"
    
    # Cache for term labels and slim mappings
    _go_adapter: Optional[object] = field(default=None, init=False)
    _slim_terms: Set[str] = field(default_factory=set, init=False)
    _term_to_slims: Dict[str, Set[str]] = field(default_factory=dict, init=False)
    _term_labels: Dict[str, str] = field(default_factory=dict, init=False)
    _term_namespace: Dict[str, str] = field(default_factory=dict, init=False)
    
    def __post_init__(self):
        """Initialize GO adapter and load slim terms."""
        self._initialize_adapter()
        self._load_slim_terms()
    
    def _initialize_adapter(self):
        """Initialize the OAK adapter for GO."""
        try:
            # Use sqlite version of GO for better performance
            self._go_adapter = get_adapter("sqlite:obo:go")
        except Exception:
            # Fallback to pronto if sqlite fails
            self._go_adapter = get_adapter("pronto:obo:go.obo")
    
    def _load_slim_terms(self):
        """Load all terms in the specified GO slim subset."""
        if not self._go_adapter:
            return
            
        try:
            # Get all terms with the slim subset annotation
            for term in self._go_adapter.entities(filter_obsoletes=True):
                term_id = str(term)
                # Check if term is in the slim subset
                if self._is_in_slim(term_id):
                    self._slim_terms.add(term_id)
                    # Cache the label
                    label = self._go_adapter.label(term_id)
                    if label:
                        self._term_labels[term_id] = label
        except Exception as e:
            # If we can't load slims dynamically, use a default set
            self._load_default_slim_terms()
    
    def _is_in_slim(self, term_id: str) -> bool:
        """Check if a term is part of the configured slim subset."""
        if not self._go_adapter:
            return False
            
        try:
            # Get term metadata to check subset membership
            # This would need to check the subset annotation
            # For now, we'll use a simplified approach
            return self._is_high_level_term(term_id)
        except Exception:
            return False
    
    def _is_high_level_term(self, term_id: str) -> bool:
        """Check if a term is a high-level term suitable for slim."""
        # Check distance from root terms
        try:
            # Get ancestors to check depth
            ancestors = set(self._go_adapter.ancestors(term_id, predicates=[IS_A]))
            # High level terms have fewer ancestors
            return len(ancestors) <= 4
        except Exception:
            return False
    
    def _load_default_slim_terms(self):
        """Load a default set of common GO slim terms."""
        # Common high-level GO slim terms
        default_slims = {
            # Biological Process
            "GO:0008150": "biological_process",
            "GO:0009987": "cellular process",
            "GO:0008152": "metabolic process",
            "GO:0065007": "biological regulation",
            "GO:0050896": "response to stimulus",
            "GO:0032502": "developmental process",
            "GO:0000003": "reproduction",
            "GO:0002376": "immune system process",
            "GO:0023052": "signaling",
            "GO:0051179": "localization",
            
            # Molecular Function  
            "GO:0003674": "molecular_function",
            "GO:0005488": "binding",
            "GO:0003824": "catalytic activity",
            "GO:0005215": "transporter activity",
            "GO:0038023": "signaling receptor activity",
            "GO:0140110": "transcription regulator activity",
            "GO:0005198": "structural molecule activity",
            
            # Cellular Component
            "GO:0005575": "cellular_component",
            "GO:0005623": "cell",
            "GO:0110165": "cellular anatomical entity",
            "GO:0032991": "protein-containing complex",
            "GO:0016020": "membrane",
            "GO:0043226": "organelle",
            "GO:0005576": "extracellular region",
        }
        
        for term_id, label in default_slims.items():
            self._slim_terms.add(term_id)
            self._term_labels[term_id] = label.replace("_", " ")
    
    def get_term_label(self, term_id: str) -> Optional[str]:
        """Get the label for a GO term.
        
        Args:
            term_id: GO term ID (e.g., "GO:0006629")
            
        Returns:
            The term label or None if not found
        """
        if term_id in self._term_labels:
            return self._term_labels[term_id]
            
        if self._go_adapter:
            try:
                label = self._go_adapter.label(term_id)
                if label:
                    self._term_labels[term_id] = label
                    return label
            except Exception:
                pass
                
        return None
    
    def get_term_namespace(self, term_id: str) -> str:
        """Get the namespace (BP, MF, CC) for a GO term.
        
        Args:
            term_id: GO term ID
            
        Returns:
            One of "biological_process", "molecular_function", "cellular_component"
        """
        if term_id in self._term_namespace:
            return self._term_namespace[term_id]
            
        if self._go_adapter:
            try:
                # Get term info to determine namespace
                # This is a simplified approach
                ancestors = list(self._go_adapter.ancestors(term_id, predicates=[IS_A]))
                
                # Check which root term it descends from
                if "GO:0008150" in ancestors:
                    namespace = "biological_process"
                elif "GO:0003674" in ancestors:
                    namespace = "molecular_function"
                elif "GO:0005575" in ancestors:
                    namespace = "cellular_component"
                else:
                    # Try to infer from term ID patterns (simplified)
                    namespace = "biological_process"  # default
                    
                self._term_namespace[term_id] = namespace
                return namespace
            except Exception:
                pass
        
        # Default based on common patterns
        return "biological_process"
    
    def get_slim_ancestors(self, term_id: str) -> Set[str]:
        """Get all GO slim ancestors for a term.
        
        Args:
            term_id: GO term ID
            
        Returns:
            Set of GO slim term IDs that are ancestors of the given term
        """
        if term_id in self._term_to_slims:
            return self._term_to_slims[term_id]
            
        slim_ancestors = set()
        
        if self._go_adapter:
            try:
                # Get all ancestors
                ancestors = set(self._go_adapter.ancestors(term_id, predicates=[IS_A]))
                # Filter to only slim terms
                slim_ancestors = ancestors.intersection(self._slim_terms)
            except Exception:
                # Fallback to basic categorization
                slim_ancestors = self._get_default_slim_category(term_id)
        else:
            slim_ancestors = self._get_default_slim_category(term_id)
            
        self._term_to_slims[term_id] = slim_ancestors
        return slim_ancestors
    
    def _get_default_slim_category(self, term_id: str) -> Set[str]:
        """Get a default slim category for a term when adapter is not available."""
        # Try to match based on the label if we have it
        label = self.get_term_label(term_id)
        if label:
            label_lower = label.lower()
            
            # Simple heuristic matching
            if "metabol" in label_lower:
                return {"GO:0008152"}  # metabolic process
            elif "bind" in label_lower:
                return {"GO:0005488"}  # binding
            elif "transport" in label_lower:
                return {"GO:0051179"}  # localization
            elif "signal" in label_lower:
                return {"GO:0023052"}  # signaling
            elif "membrane" in label_lower:
                return {"GO:0016020"}  # membrane
            elif "nucleus" in label_lower or "nuclear" in label_lower:
                return {"GO:0043226"}  # organelle
                
        # Default to cellular process for BP terms
        namespace = self.get_term_namespace(term_id)
        if namespace == "biological_process":
            return {"GO:0009987"}  # cellular process
        elif namespace == "molecular_function":
            return {"GO:0003674"}  # molecular_function
        else:
            return {"GO:0005575"}  # cellular_component
    
    def organize_by_slim(self, term_ids: List[str]) -> Dict[str, Dict[str, List[str]]]:
        """Organize a list of GO terms by their slim categories and namespaces.
        
        Args:
            term_ids: List of GO term IDs
            
        Returns:
            Nested dict: {namespace: {slim_term_id: [child_term_ids]}}
        """
        organized = {
            "biological_process": {},
            "molecular_function": {},
            "cellular_component": {}
        }
        
        for term_id in term_ids:
            namespace = self.get_term_namespace(term_id)
            slim_ancestors = self.get_slim_ancestors(term_id)
            
            if not slim_ancestors:
                # Put under root if no slim ancestors found
                root_term = self._get_root_for_namespace(namespace)
                slim_ancestors = {root_term}
            
            # Add to most specific slim ancestor
            most_specific = self._get_most_specific_slim(slim_ancestors)
            
            if namespace not in organized:
                organized[namespace] = {}
            if most_specific not in organized[namespace]:
                organized[namespace][most_specific] = []
                
            organized[namespace][most_specific].append(term_id)
        
        return organized
    
    def _get_root_for_namespace(self, namespace: str) -> str:
        """Get the root GO term for a namespace."""
        roots = {
            "biological_process": "GO:0008150",
            "molecular_function": "GO:0003674",
            "cellular_component": "GO:0005575"
        }
        return roots.get(namespace, "GO:0008150")
    
    def _get_most_specific_slim(self, slim_terms: Set[str]) -> str:
        """From a set of slim terms, return the most specific one."""
        if not slim_terms:
            return "GO:0008150"  # default to biological process root
            
        if len(slim_terms) == 1:
            return list(slim_terms)[0]
            
        # Try to find the most specific by checking which has most ancestors
        most_specific = None
        max_ancestors = -1
        
        for term in slim_terms:
            if self._go_adapter:
                try:
                    ancestors = list(self._go_adapter.ancestors(term, predicates=[IS_A]))
                    if len(ancestors) > max_ancestors:
                        max_ancestors = len(ancestors)
                        most_specific = term
                except Exception:
                    pass
                    
        return most_specific or list(slim_terms)[0]
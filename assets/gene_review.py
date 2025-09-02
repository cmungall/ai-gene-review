# Auto generated from gene_review.yaml by pythongen.py version: 0.0.1
# Generation date: 2025-08-31T14:45:56
# Schema: gene_curation
#
# id: https://w3id.org/ai4curation/gene_review
# description: Schema for gene curation Top level entity is a GeneReview, which is about a single gene (and its equivalent swiss-prot entry). It contains a high level summary of the gene, plus a review of all existing annotations. It also contains a list of core functions, which are GO-CAM-like annotons describing the core evolved functions of the gene.
# license: https://creativecommons.org/publicdomain/zero/1.0/

from dataclasses import dataclass
from typing import (
    Any,
    ClassVar,
    Optional,
    Union
)

from jsonasobj2 import (
    as_dict
)
from linkml_runtime.linkml_model.meta import (
    EnumDefinition,
    PermissibleValue
)
from linkml_runtime.utils.curienamespace import CurieNamespace
from linkml_runtime.utils.enumerations import EnumDefinitionImpl
from linkml_runtime.utils.metamodelcore import (
    empty_dict,
    empty_list
)
from linkml_runtime.utils.slot import Slot
from linkml_runtime.utils.yamlutils import (
    YAMLRoot,
    extended_str
)
from rdflib import (
    URIRef
)

from linkml_runtime.utils.metamodelcore import Bool

metamodel_version = "1.7.0"
version = None

# Namespaces
DCAT = CurieNamespace('dcat', 'http://www.w3.org/ns/dcat#')
DCTERMS = CurieNamespace('dcterms', 'http://purl.org/dc/terms/')
GENE_REVIEW = CurieNamespace('gene_review', 'https://w3id.org/ai4curation/gene_review/')
LINKML = CurieNamespace('linkml', 'https://w3id.org/linkml/')
OWL = CurieNamespace('owl', 'http://www.w3.org/2002/07/owl#')
RDF = CurieNamespace('rdf', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#')
RDFS = CurieNamespace('rdfs', 'http://www.w3.org/2000/01/rdf-schema#')
SKOS = CurieNamespace('skos', 'http://www.w3.org/2004/02/skos/core#')
XSD = CurieNamespace('xsd', 'http://www.w3.org/2001/XMLSchema#')
DEFAULT_ = GENE_REVIEW


# Types

# Class references
class GeneReviewId(extended_str):
    pass


class TermId(extended_str):
    pass


class ReferenceId(extended_str):
    pass


@dataclass(repr=False)
class GeneReview(YAMLRoot):
    """
    Complete review for a gene
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["GeneReview"]
    class_class_curie: ClassVar[str] = "gene_review:GeneReview"
    class_name: ClassVar[str] = "GeneReview"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.GeneReview

    id: Union[str, GeneReviewId] = None
    gene_symbol: str = None
    aliases: Optional[Union[str, list[str]]] = empty_list()
    description: Optional[str] = None
    taxon: Optional[Union[dict, "Term"]] = None
    references: Optional[Union[dict[Union[str, ReferenceId], Union[dict, "Reference"]], list[Union[dict, "Reference"]]]] = empty_dict()
    existing_annotations: Optional[Union[Union[dict, "ExistingAnnotation"], list[Union[dict, "ExistingAnnotation"]]]] = empty_list()
    core_functions: Optional[Union[Union[dict, "CoreFunction"], list[Union[dict, "CoreFunction"]]]] = empty_list()
    proposed_new_terms: Optional[Union[Union[dict, "ProposedOntologyTerm"], list[Union[dict, "ProposedOntologyTerm"]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, GeneReviewId):
            self.id = GeneReviewId(self.id)

        if self._is_empty(self.gene_symbol):
            self.MissingRequiredField("gene_symbol")
        if not isinstance(self.gene_symbol, str):
            self.gene_symbol = str(self.gene_symbol)

        if not isinstance(self.aliases, list):
            self.aliases = [self.aliases] if self.aliases is not None else []
        self.aliases = [v if isinstance(v, str) else str(v) for v in self.aliases]

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.taxon is not None and not isinstance(self.taxon, Term):
            self.taxon = Term(**as_dict(self.taxon))

        self._normalize_inlined_as_list(slot_name="references", slot_type=Reference, key_name="id", keyed=True)

        if not isinstance(self.existing_annotations, list):
            self.existing_annotations = [self.existing_annotations] if self.existing_annotations is not None else []
        self.existing_annotations = [v if isinstance(v, ExistingAnnotation) else ExistingAnnotation(**as_dict(v)) for v in self.existing_annotations]

        if not isinstance(self.core_functions, list):
            self.core_functions = [self.core_functions] if self.core_functions is not None else []
        self.core_functions = [v if isinstance(v, CoreFunction) else CoreFunction(**as_dict(v)) for v in self.core_functions]

        self._normalize_inlined_as_dict(slot_name="proposed_new_terms", slot_type=ProposedOntologyTerm, key_name="proposed_name", keyed=False)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Term(YAMLRoot):
    """
    A term in a specific ontology
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["Term"]
    class_class_curie: ClassVar[str] = "gene_review:Term"
    class_name: ClassVar[str] = "Term"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.Term

    id: Union[str, TermId] = None
    label: str = None
    description: Optional[str] = None
    ontology: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, TermId):
            self.id = TermId(self.id)

        if self._is_empty(self.label):
            self.MissingRequiredField("label")
        if not isinstance(self.label, str):
            self.label = str(self.label)

        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if self.ontology is not None and not isinstance(self.ontology, str):
            self.ontology = str(self.ontology)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Reference(YAMLRoot):
    """
    A reference is a published text that describes a finding or a method. References might be formal publications
    (where the ID is a PMID), or for methods, a GO_REF. Additionally, a reference to a local ad-hoc analysis or review
    can be made by using the `file:` prefix.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["Reference"]
    class_class_curie: ClassVar[str] = "gene_review:Reference"
    class_name: ClassVar[str] = "Reference"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.Reference

    id: Union[str, ReferenceId] = None
    title: str = None
    findings: Optional[Union[Union[dict, "Finding"], list[Union[dict, "Finding"]]]] = empty_list()
    is_invalid: Optional[Union[bool, Bool]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.id):
            self.MissingRequiredField("id")
        if not isinstance(self.id, ReferenceId):
            self.id = ReferenceId(self.id)

        if self._is_empty(self.title):
            self.MissingRequiredField("title")
        if not isinstance(self.title, str):
            self.title = str(self.title)

        if not isinstance(self.findings, list):
            self.findings = [self.findings] if self.findings is not None else []
        self.findings = [v if isinstance(v, Finding) else Finding(**as_dict(v)) for v in self.findings]

        if self.is_invalid is not None and not isinstance(self.is_invalid, Bool):
            self.is_invalid = Bool(self.is_invalid)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Finding(YAMLRoot):
    """
    A finding is a statement about a gene, which is supported by a reference. Similar to "comments" in uniprot
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["Finding"]
    class_class_curie: ClassVar[str] = "gene_review:Finding"
    class_name: ClassVar[str] = "Finding"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.Finding

    statement: Optional[str] = None
    supporting_text: Optional[str] = None
    reference_section_type: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.statement is not None and not isinstance(self.statement, str):
            self.statement = str(self.statement)

        if self.supporting_text is not None and not isinstance(self.supporting_text, str):
            self.supporting_text = str(self.supporting_text)

        if self.reference_section_type is not None and not isinstance(self.reference_section_type, str):
            self.reference_section_type = str(self.reference_section_type)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class SupportingTextInReference(YAMLRoot):
    """
    A supporting text in a reference.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["SupportingTextInReference"]
    class_class_curie: ClassVar[str] = "gene_review:SupportingTextInReference"
    class_name: ClassVar[str] = "SupportingTextInReference"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.SupportingTextInReference

    reference_id: Optional[Union[str, ReferenceId]] = None
    supporting_text: Optional[str] = None
    reference_section_type: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.reference_id is not None and not isinstance(self.reference_id, ReferenceId):
            self.reference_id = ReferenceId(self.reference_id)

        if self.supporting_text is not None and not isinstance(self.supporting_text, str):
            self.supporting_text = str(self.supporting_text)

        if self.reference_section_type is not None and not isinstance(self.reference_section_type, str):
            self.reference_section_type = str(self.reference_section_type)

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ExistingAnnotation(YAMLRoot):
    """
    An existing annotation from the GO database, plus a review of the annotation.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["ExistingAnnotation"]
    class_class_curie: ClassVar[str] = "gene_review:ExistingAnnotation"
    class_name: ClassVar[str] = "ExistingAnnotation"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.ExistingAnnotation

    term: Optional[Union[dict, Term]] = None
    extensions: Optional[Union[Union[dict, "AnnotationExtension"], list[Union[dict, "AnnotationExtension"]]]] = empty_list()
    negated: Optional[Union[bool, Bool]] = None
    evidence_type: Optional[str] = None
    original_reference_id: Optional[Union[str, ReferenceId]] = None
    supporting_entities: Optional[Union[str, list[str]]] = empty_list()
    review: Optional[Union[dict, "Review"]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.term is not None and not isinstance(self.term, Term):
            self.term = Term(**as_dict(self.term))

        if not isinstance(self.extensions, list):
            self.extensions = [self.extensions] if self.extensions is not None else []
        self.extensions = [v if isinstance(v, AnnotationExtension) else AnnotationExtension(**as_dict(v)) for v in self.extensions]

        if self.negated is not None and not isinstance(self.negated, Bool):
            self.negated = Bool(self.negated)

        if self.evidence_type is not None and not isinstance(self.evidence_type, str):
            self.evidence_type = str(self.evidence_type)

        if self.original_reference_id is not None and not isinstance(self.original_reference_id, ReferenceId):
            self.original_reference_id = ReferenceId(self.original_reference_id)

        if not isinstance(self.supporting_entities, list):
            self.supporting_entities = [self.supporting_entities] if self.supporting_entities is not None else []
        self.supporting_entities = [v if isinstance(v, str) else str(v) for v in self.supporting_entities]

        if self.review is not None and not isinstance(self.review, Review):
            self.review = Review(**as_dict(self.review))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class Review(YAMLRoot):
    """
    A review of an existing annotation.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["Review"]
    class_class_curie: ClassVar[str] = "gene_review:Review"
    class_name: ClassVar[str] = "Review"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.Review

    summary: Optional[str] = None
    action: Optional[Union[str, "ActionEnum"]] = None
    reason: Optional[str] = None
    proposed_replacement_terms: Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]] = empty_dict()
    additional_reference_ids: Optional[Union[Union[str, ReferenceId], list[Union[str, ReferenceId]]]] = empty_list()
    supported_by: Optional[Union[Union[dict, SupportingTextInReference], list[Union[dict, SupportingTextInReference]]]] = empty_list()

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.summary is not None and not isinstance(self.summary, str):
            self.summary = str(self.summary)

        if self.action is not None and not isinstance(self.action, ActionEnum):
            self.action = ActionEnum(self.action)

        if self.reason is not None and not isinstance(self.reason, str):
            self.reason = str(self.reason)

        self._normalize_inlined_as_list(slot_name="proposed_replacement_terms", slot_type=Term, key_name="id", keyed=True)

        if not isinstance(self.additional_reference_ids, list):
            self.additional_reference_ids = [self.additional_reference_ids] if self.additional_reference_ids is not None else []
        self.additional_reference_ids = [v if isinstance(v, ReferenceId) else ReferenceId(v) for v in self.additional_reference_ids]

        if not isinstance(self.supported_by, list):
            self.supported_by = [self.supported_by] if self.supported_by is not None else []
        self.supported_by = [v if isinstance(v, SupportingTextInReference) else SupportingTextInReference(**as_dict(v)) for v in self.supported_by]

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class CoreFunction(YAMLRoot):
    """
    A core function is a GO-CAM-like annotation of the core evolved functions of a gene. This is a synthesis of the
    reviewed core annotations, brought together into a unified GO-CAM-like representation.
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["CoreFunction"]
    class_class_curie: ClassVar[str] = "gene_review:CoreFunction"
    class_name: ClassVar[str] = "CoreFunction"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.CoreFunction

    description: Optional[str] = None
    supported_by: Optional[Union[Union[dict, SupportingTextInReference], list[Union[dict, SupportingTextInReference]]]] = empty_list()
    molecular_function: Optional[Union[dict, Term]] = None
    directly_involved_in: Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]] = empty_dict()
    locations: Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]] = empty_dict()
    anatomical_locations: Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]] = empty_dict()
    substrates: Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]] = empty_dict()
    in_complex: Optional[Union[dict, Term]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.description is not None and not isinstance(self.description, str):
            self.description = str(self.description)

        if not isinstance(self.supported_by, list):
            self.supported_by = [self.supported_by] if self.supported_by is not None else []
        self.supported_by = [v if isinstance(v, SupportingTextInReference) else SupportingTextInReference(**as_dict(v)) for v in self.supported_by]

        if self.molecular_function is not None and not isinstance(self.molecular_function, Term):
            self.molecular_function = Term(**as_dict(self.molecular_function))

        self._normalize_inlined_as_list(slot_name="directly_involved_in", slot_type=Term, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="locations", slot_type=Term, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="anatomical_locations", slot_type=Term, key_name="id", keyed=True)

        self._normalize_inlined_as_list(slot_name="substrates", slot_type=Term, key_name="id", keyed=True)

        if self.in_complex is not None and not isinstance(self.in_complex, Term):
            self.in_complex = Term(**as_dict(self.in_complex))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class AnnotationExtension(YAMLRoot):
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["AnnotationExtension"]
    class_class_curie: ClassVar[str] = "gene_review:AnnotationExtension"
    class_name: ClassVar[str] = "AnnotationExtension"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.AnnotationExtension

    predicate: Optional[str] = None
    term: Optional[Union[dict, Term]] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self.predicate is not None and not isinstance(self.predicate, str):
            self.predicate = str(self.predicate)

        if self.term is not None and not isinstance(self.term, Term):
            self.term = Term(**as_dict(self.term))

        super().__post_init__(**kwargs)


@dataclass(repr=False)
class ProposedOntologyTerm(YAMLRoot):
    """
    A proposed new ontology term that should exist but doesn't currently
    """
    _inherited_slots: ClassVar[list[str]] = []

    class_class_uri: ClassVar[URIRef] = GENE_REVIEW["ProposedOntologyTerm"]
    class_class_curie: ClassVar[str] = "gene_review:ProposedOntologyTerm"
    class_name: ClassVar[str] = "ProposedOntologyTerm"
    class_model_uri: ClassVar[URIRef] = GENE_REVIEW.ProposedOntologyTerm

    proposed_name: str = None
    proposed_definition: str = None
    justification: Optional[str] = None

    def __post_init__(self, *_: str, **kwargs: Any):
        if self._is_empty(self.proposed_name):
            self.MissingRequiredField("proposed_name")
        if not isinstance(self.proposed_name, str):
            self.proposed_name = str(self.proposed_name)

        if self._is_empty(self.proposed_definition):
            self.MissingRequiredField("proposed_definition")
        if not isinstance(self.proposed_definition, str):
            self.proposed_definition = str(self.proposed_definition)

        if self.justification is not None and not isinstance(self.justification, str):
            self.justification = str(self.justification)

        super().__post_init__(**kwargs)


# Enumerations
class EvidenceType(EnumDefinitionImpl):

    IDA = PermissibleValue(text="IDA")
    IBA = PermissibleValue(text="IBA")
    ISS = PermissibleValue(text="ISS")
    TAS = PermissibleValue(text="TAS")
    IEP = PermissibleValue(text="IEP")
    IGC = PermissibleValue(text="IGC")

    _defn = EnumDefinition(
        name="EvidenceType",
    )

class ActionEnum(EnumDefinitionImpl):

    ACCEPT = PermissibleValue(
        text="ACCEPT",
        description="""Accept the existing annotation as-is, no modifications, and retain as representing the core function of the gene""")
    KEEP_AS_NON_CORE = PermissibleValue(
        text="KEEP_AS_NON_CORE",
        description="""Keep the existing annotation as-is, but mark it as non-core. For pleiotropic genes, this may be the developmental processes, or other processes that are not the core function of the gene.""")
    REMOVE = PermissibleValue(
        text="REMOVE",
        description="Remove the existing annotation, as it is unlikely to be correct based on combined evidence")
    MODIFY = PermissibleValue(
        text="MODIFY",
        description="""The essence of the annotation is sound, but there are better terms to use (use in combination with proposed_replacement_terms). if the term is too general, then MODIFY should be used, with a proposed replacement term for the correct specific function. sometimes terms can also be overly specific and contorted, so in some cases you might want to generalize""")
    MARK_AS_OVER_ANNOTATED = PermissibleValue(
        text="MARK_AS_OVER_ANNOTATED",
        description="The term is not entirely wrong, but likely represents an over-annotation of the gene")
    UNDECIDED = PermissibleValue(
        text="UNDECIDED",
        description="""The annotation is not clear, and the reviewer is not sure what to do with it. ALWAYS USE THIS IF YOU ARE UNABLE TO ACCESS RELEVANT PUBLICATIONS""")
    PENDING = PermissibleValue(
        text="PENDING",
        description="The review entry is a stub, and the review has not been completed yet.")

    _defn = EnumDefinition(
        name="ActionEnum",
    )

class GOTermEnum(EnumDefinitionImpl):
    """
    A term in the GO ontology
    """
    _defn = EnumDefinition(
        name="GOTermEnum",
        description="A term in the GO ontology",
    )

class GOMolecularActivityEnum(EnumDefinitionImpl):
    """
    A molecular activity term in the GO ontology
    """
    _defn = EnumDefinition(
        name="GOMolecularActivityEnum",
        description="A molecular activity term in the GO ontology",
    )

class GOBiologicalProcessEnum(EnumDefinitionImpl):
    """
    A biological process term in the GO ontology
    """
    _defn = EnumDefinition(
        name="GOBiologicalProcessEnum",
        description="A biological process term in the GO ontology",
    )

class GOCellularLocationEnum(EnumDefinitionImpl):
    """
    A cellular location term in the GO ontology (excludes protein-containing complexes)
    """
    _defn = EnumDefinition(
        name="GOCellularLocationEnum",
        description="A cellular location term in the GO ontology (excludes protein-containing complexes)",
    )

class GOProteinContainingComplexEnum(EnumDefinitionImpl):
    """
    A protein-containing complex term in the GO ontology
    """
    _defn = EnumDefinition(
        name="GOProteinContainingComplexEnum",
        description="A protein-containing complex term in the GO ontology",
    )

# Slots
class slots:
    pass

slots.id = Slot(uri=GENE_REVIEW.id, name="id", curie=GENE_REVIEW.curie('id'),
                   model_uri=GENE_REVIEW.id, domain=None, range=URIRef)

slots.label = Slot(uri=GENE_REVIEW.label, name="label", curie=GENE_REVIEW.curie('label'),
                   model_uri=GENE_REVIEW.label, domain=None, range=str)

slots.gene_symbol = Slot(uri=GENE_REVIEW.gene_symbol, name="gene_symbol", curie=GENE_REVIEW.curie('gene_symbol'),
                   model_uri=GENE_REVIEW.gene_symbol, domain=None, range=str)

slots.title = Slot(uri=DCTERMS.title, name="title", curie=DCTERMS.curie('title'),
                   model_uri=GENE_REVIEW.title, domain=None, range=str)

slots.aliases = Slot(uri=GENE_REVIEW.aliases, name="aliases", curie=GENE_REVIEW.curie('aliases'),
                   model_uri=GENE_REVIEW.aliases, domain=None, range=Optional[Union[str, list[str]]])

slots.description = Slot(uri=GENE_REVIEW.description, name="description", curie=GENE_REVIEW.curie('description'),
                   model_uri=GENE_REVIEW.description, domain=None, range=Optional[str])

slots.statement = Slot(uri=GENE_REVIEW.statement, name="statement", curie=GENE_REVIEW.curie('statement'),
                   model_uri=GENE_REVIEW.statement, domain=None, range=Optional[str])

slots.references = Slot(uri=GENE_REVIEW.references, name="references", curie=GENE_REVIEW.curie('references'),
                   model_uri=GENE_REVIEW.references, domain=None, range=Optional[Union[dict[Union[str, ReferenceId], Union[dict, Reference]], list[Union[dict, Reference]]]])

slots.findings = Slot(uri=GENE_REVIEW.findings, name="findings", curie=GENE_REVIEW.curie('findings'),
                   model_uri=GENE_REVIEW.findings, domain=None, range=Optional[Union[Union[dict, Finding], list[Union[dict, Finding]]]])

slots.supported_by = Slot(uri=GENE_REVIEW.supported_by, name="supported_by", curie=GENE_REVIEW.curie('supported_by'),
                   model_uri=GENE_REVIEW.supported_by, domain=None, range=Optional[Union[Union[dict, SupportingTextInReference], list[Union[dict, SupportingTextInReference]]]])

slots.summary = Slot(uri=GENE_REVIEW.summary, name="summary", curie=GENE_REVIEW.curie('summary'),
                   model_uri=GENE_REVIEW.summary, domain=None, range=Optional[str])

slots.supporting_text = Slot(uri=GENE_REVIEW.supporting_text, name="supporting_text", curie=GENE_REVIEW.curie('supporting_text'),
                   model_uri=GENE_REVIEW.supporting_text, domain=None, range=Optional[str])

slots.evidence_type = Slot(uri=GENE_REVIEW.evidence_type, name="evidence_type", curie=GENE_REVIEW.curie('evidence_type'),
                   model_uri=GENE_REVIEW.evidence_type, domain=None, range=Optional[str])

slots.term = Slot(uri=GENE_REVIEW.term, name="term", curie=GENE_REVIEW.curie('term'),
                   model_uri=GENE_REVIEW.term, domain=None, range=Optional[Union[dict, Term]])

slots.predicate = Slot(uri=GENE_REVIEW.predicate, name="predicate", curie=GENE_REVIEW.curie('predicate'),
                   model_uri=GENE_REVIEW.predicate, domain=None, range=Optional[str])

slots.taxon = Slot(uri=GENE_REVIEW.taxon, name="taxon", curie=GENE_REVIEW.curie('taxon'),
                   model_uri=GENE_REVIEW.taxon, domain=None, range=Optional[Union[dict, Term]])

slots.existing_annotations = Slot(uri=GENE_REVIEW.existing_annotations, name="existing_annotations", curie=GENE_REVIEW.curie('existing_annotations'),
                   model_uri=GENE_REVIEW.existing_annotations, domain=None, range=Optional[Union[Union[dict, ExistingAnnotation], list[Union[dict, ExistingAnnotation]]]])

slots.core_functions = Slot(uri=GENE_REVIEW.core_functions, name="core_functions", curie=GENE_REVIEW.curie('core_functions'),
                   model_uri=GENE_REVIEW.core_functions, domain=None, range=Optional[Union[Union[dict, CoreFunction], list[Union[dict, CoreFunction]]]])

slots.action = Slot(uri=GENE_REVIEW.action, name="action", curie=GENE_REVIEW.curie('action'),
                   model_uri=GENE_REVIEW.action, domain=None, range=Optional[Union[str, "ActionEnum"]])

slots.reason = Slot(uri=GENE_REVIEW.reason, name="reason", curie=GENE_REVIEW.curie('reason'),
                   model_uri=GENE_REVIEW.reason, domain=None, range=Optional[str])

slots.proposed_replacement_terms = Slot(uri=GENE_REVIEW.proposed_replacement_terms, name="proposed_replacement_terms", curie=GENE_REVIEW.curie('proposed_replacement_terms'),
                   model_uri=GENE_REVIEW.proposed_replacement_terms, domain=None, range=Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]])

slots.extensions = Slot(uri=GENE_REVIEW.extensions, name="extensions", curie=GENE_REVIEW.curie('extensions'),
                   model_uri=GENE_REVIEW.extensions, domain=None, range=Optional[Union[Union[dict, AnnotationExtension], list[Union[dict, AnnotationExtension]]]])

slots.negated = Slot(uri=GENE_REVIEW.negated, name="negated", curie=GENE_REVIEW.curie('negated'),
                   model_uri=GENE_REVIEW.negated, domain=None, range=Optional[Union[bool, Bool]])

slots.reference_id = Slot(uri=GENE_REVIEW.reference_id, name="reference_id", curie=GENE_REVIEW.curie('reference_id'),
                   model_uri=GENE_REVIEW.reference_id, domain=None, range=Optional[Union[str, ReferenceId]])

slots.original_reference_id = Slot(uri=GENE_REVIEW.original_reference_id, name="original_reference_id", curie=GENE_REVIEW.curie('original_reference_id'),
                   model_uri=GENE_REVIEW.original_reference_id, domain=None, range=Optional[Union[str, ReferenceId]])

slots.review = Slot(uri=GENE_REVIEW.review, name="review", curie=GENE_REVIEW.curie('review'),
                   model_uri=GENE_REVIEW.review, domain=None, range=Optional[Union[dict, Review]])

slots.ontology = Slot(uri=GENE_REVIEW.ontology, name="ontology", curie=GENE_REVIEW.curie('ontology'),
                   model_uri=GENE_REVIEW.ontology, domain=None, range=Optional[str])

slots.supporting_entities = Slot(uri=GENE_REVIEW.supporting_entities, name="supporting_entities", curie=GENE_REVIEW.curie('supporting_entities'),
                   model_uri=GENE_REVIEW.supporting_entities, domain=None, range=Optional[Union[str, list[str]]])

slots.additional_reference_ids = Slot(uri=GENE_REVIEW.additional_reference_ids, name="additional_reference_ids", curie=GENE_REVIEW.curie('additional_reference_ids'),
                   model_uri=GENE_REVIEW.additional_reference_ids, domain=None, range=Optional[Union[Union[str, ReferenceId], list[Union[str, ReferenceId]]]])

slots.is_invalid = Slot(uri=GENE_REVIEW.is_invalid, name="is_invalid", curie=GENE_REVIEW.curie('is_invalid'),
                   model_uri=GENE_REVIEW.is_invalid, domain=None, range=Optional[Union[bool, Bool]])

slots.reference_section_type = Slot(uri=GENE_REVIEW.reference_section_type, name="reference_section_type", curie=GENE_REVIEW.curie('reference_section_type'),
                   model_uri=GENE_REVIEW.reference_section_type, domain=None, range=Optional[str])

slots.proposed_new_terms = Slot(uri=GENE_REVIEW.proposed_new_terms, name="proposed_new_terms", curie=GENE_REVIEW.curie('proposed_new_terms'),
                   model_uri=GENE_REVIEW.proposed_new_terms, domain=None, range=Optional[Union[Union[dict, ProposedOntologyTerm], list[Union[dict, ProposedOntologyTerm]]]])

slots.coreFunction__description = Slot(uri=GENE_REVIEW.description, name="coreFunction__description", curie=GENE_REVIEW.curie('description'),
                   model_uri=GENE_REVIEW.coreFunction__description, domain=None, range=Optional[str])

slots.coreFunction__supported_by = Slot(uri=GENE_REVIEW.supported_by, name="coreFunction__supported_by", curie=GENE_REVIEW.curie('supported_by'),
                   model_uri=GENE_REVIEW.coreFunction__supported_by, domain=None, range=Optional[Union[Union[dict, SupportingTextInReference], list[Union[dict, SupportingTextInReference]]]])

slots.coreFunction__molecular_function = Slot(uri=GENE_REVIEW.molecular_function, name="coreFunction__molecular_function", curie=GENE_REVIEW.curie('molecular_function'),
                   model_uri=GENE_REVIEW.coreFunction__molecular_function, domain=None, range=Optional[Union[dict, Term]])

slots.coreFunction__directly_involved_in = Slot(uri=GENE_REVIEW.directly_involved_in, name="coreFunction__directly_involved_in", curie=GENE_REVIEW.curie('directly_involved_in'),
                   model_uri=GENE_REVIEW.coreFunction__directly_involved_in, domain=None, range=Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]])

slots.coreFunction__locations = Slot(uri=GENE_REVIEW.locations, name="coreFunction__locations", curie=GENE_REVIEW.curie('locations'),
                   model_uri=GENE_REVIEW.coreFunction__locations, domain=None, range=Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]])

slots.coreFunction__anatomical_locations = Slot(uri=GENE_REVIEW.anatomical_locations, name="coreFunction__anatomical_locations", curie=GENE_REVIEW.curie('anatomical_locations'),
                   model_uri=GENE_REVIEW.coreFunction__anatomical_locations, domain=None, range=Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]])

slots.coreFunction__substrates = Slot(uri=GENE_REVIEW.substrates, name="coreFunction__substrates", curie=GENE_REVIEW.curie('substrates'),
                   model_uri=GENE_REVIEW.coreFunction__substrates, domain=None, range=Optional[Union[dict[Union[str, TermId], Union[dict, Term]], list[Union[dict, Term]]]])

slots.coreFunction__in_complex = Slot(uri=GENE_REVIEW.in_complex, name="coreFunction__in_complex", curie=GENE_REVIEW.curie('in_complex'),
                   model_uri=GENE_REVIEW.coreFunction__in_complex, domain=None, range=Optional[Union[dict, Term]])

slots.proposedOntologyTerm__proposed_name = Slot(uri=GENE_REVIEW.proposed_name, name="proposedOntologyTerm__proposed_name", curie=GENE_REVIEW.curie('proposed_name'),
                   model_uri=GENE_REVIEW.proposedOntologyTerm__proposed_name, domain=None, range=str)

slots.proposedOntologyTerm__proposed_definition = Slot(uri=GENE_REVIEW.proposed_definition, name="proposedOntologyTerm__proposed_definition", curie=GENE_REVIEW.curie('proposed_definition'),
                   model_uri=GENE_REVIEW.proposedOntologyTerm__proposed_definition, domain=None, range=str)

slots.proposedOntologyTerm__justification = Slot(uri=GENE_REVIEW.justification, name="proposedOntologyTerm__justification", curie=GENE_REVIEW.curie('justification'),
                   model_uri=GENE_REVIEW.proposedOntologyTerm__justification, domain=None, range=Optional[str])

slots.ExistingAnnotation_term = Slot(uri=GENE_REVIEW.term, name="ExistingAnnotation_term", curie=GENE_REVIEW.curie('term'),
                   model_uri=GENE_REVIEW.ExistingAnnotation_term, domain=ExistingAnnotation, range=Optional[Union[dict, Term]])

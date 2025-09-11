
# Enum: ActionEnum



URI: [gene_review:ActionEnum](https://w3id.org/ai4curation/gene_review/ActionEnum)


## Permissible Values

| Text | Description | Meaning | Other Information |
| :--- | :---: | :---: | ---: |
| ACCEPT | Accept the existing annotation as-is, no modifications, and retain as representing the core function of the gene |  |  |
| KEEP_AS_NON_CORE | Keep the existing annotation as-is, but mark it as non-core. For pleiotropic genes, this may be the developmental processes, or other processes that are not the core function of the gene. |  |  |
| REMOVE | Remove the existing annotation, as it is unlikely to be correct based on combined evidence |  |  |
| MODIFY | The essence of the annotation is sound, but there are better terms to use (use in combination with proposed_replacement_terms). if the term is too general, then MODIFY should be used, with a proposed replacement term for the correct specific function. sometimes terms can also be overly specific and contorted, so in some cases you might want to generalize |  |  |
| MARK_AS_OVER_ANNOTATED | The term is not entirely wrong, but likely represents an over-annotation of the gene |  |  |
| UNDECIDED | The annotation is not clear, and the reviewer is not sure what to do with it. ALWAYS USE THIS IF YOU ARE UNABLE TO ACCESS RELEVANT PUBLICATIONS |  |  |
| PENDING | The review entry is a stub, and the review has not been completed yet. |  |  |
| NEW | This is a proposed annotation, not one that exists in the existing GO annotations. Use this to propose a new annotation not covered by the existing GO annotations. Use this conservatively, do not over-annotate, especially for biological process. Do not use for indirect or pleiotropic effects. Be sure you have good evidence, this can be from multiple sources. |  |  |

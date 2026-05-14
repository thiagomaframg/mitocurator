# Modelo de decisão da pipeline (evidence-driven)

## Princípio
A curadoria deve seguir um fluxo rastreável:

**Problema observado → Evidência → Diagnóstico → Proposta de correção → Validação → Relatório**

Este modelo descreve como os módulos atuais do MitoCurator implementam (total ou parcialmente) essa lógica.

---

## 1) Problema observado
Origem do problema no estado atual:
- `gene_qc.tsv` e `problematic_features.tsv`:
  - `internal_stop_count > 0`
  - `multiple_of_three = no`
  - combinação dos dois sinais
- `intergenic_regions.tsv`:
  - lacunas longas/AT-rich sugestivas de regiões biológicas relevantes
- `expected_gene_set.tsv` (refinement):
  - genes `MISSING` ou `DUPLICATED`

## 2) Evidência
Evidências computadas automaticamente:
- Coordenadas, strands, comprimentos nt/aa.
- Posições de stop interno.
- %AT de lacunas intergênicas.
- Ranking de candidatos de ORF para genes ausentes.
- Candidatos de ajuste local para CDS problemáticas.

Evidências ainda não cobertas no snapshot:
- suporte de leitura (BAM/coverage por gene)
- consenso direcionado por região
- comparação multi-assembly por pools

## 3) Diagnóstico
Regras diagnósticas explícitas já codificadas:
- `CHECK_INTERNAL_STOP`
- `CHECK_LENGTH_NOT_MULTIPLE_OF_THREE`
- `CHECK_FRAMESHIFT_AND_INTERNAL_STOP`
- hints de refinamento: `STRONG_CANDIDATE`, `SUGGEST_REVIEW`, `LOW_PRIORITY`, etc.

Interpretação prática:
- Problemas de frame/stop indicam risco de truncamento, indel, erro de fronteira CDS ou pseudogene.
- Ausência de gene esperado indica lacuna de anotação ou montagem incompleta/localmente incorreta.

## 4) Proposta de correção
No escopo atual, propostas são **assistidas por heurística**, não aplicadas cegamente:
- revisão manual de fronteiras CDS com base em candidatos do refinamento;
- revisão manual de regiões intergênicas para possíveis genes ausentes;
- anotação de região AT-rich putativa quando critérios de tamanho/%AT são atendidos.

## 5) Validação
Com os módulos hoje disponíveis, validação mínima recomendada:
1. Reexecutar `diagnose` após ajuste manual da anotação.
2. Confirmar redução de sinais de problema (`problematic_features.tsv`).
3. Verificar coerência global no `diagnostic_summary.md`.

Validação robusta planejada (fora do snapshot):
- confirmar candidatos com mapeamento de leituras;
- checar consenso local por extração direcionada;
- comparar assemblies/pools candidatos.

## 6) Relatório
Artefatos de relatório já gerados:
- `07_gene_qc/diagnostic_summary.md`
- tabelas TSV de diagnóstico e refinamento

Recomendação de governança curatorial:
- Cada decisão manual deve referenciar: arquivo TSV de origem, coordenadas, motivo da decisão e resultado pós-validação.

---

## Current evidence-driven execution model

The MitoCurator workflow is organized as an evidence-driven curation pipeline. The intended order is:

1. **Input mitochondrial sequence**
   - FASTA sequence recovered by tools such as MitoHiFi, PipeASM or another assembly/recovery workflow; or
   - an already annotated GenBank file.

2. **Initial annotation**
   - If the input is FASTA, MitoFinder is used to generate an initial GenBank annotation.
   - If the input is GenBank, the existing annotation is used as the starting point.

3. **Initial refinement**
   - Standardize gene names and feature labels.
   - Compare the annotation against an expected mitochondrial gene set.
   - Add or normalize non-coding regions such as AT-rich/control-region features when supported.
   - Generate refinement candidates for missing genes or problematic CDS features.

4. **Diagnosis and annotation assessment**
   - Measure annotation completeness and feature-level quality.
   - Detect missing genes, problematic CDS features, internal stops, possible boundary issues, and missing tRNAs/rRNAs.
   - Generate a curator-facing annotation assessment report.
   - Outputs include:
     - `annotation_gene_inventory.tsv`
     - `annotation_review_targets.tsv`
     - `annotation_evidence_summary.tsv`
     - `annotation_assessment_report.md`
     - `annotation_assessment_report.html`

5. **Read mapping to the mitogenome**
   - Map available sequencing reads against the current mitochondrial reference.
   - Generate sorted/indexed BAM files, coverage summaries and mapping statistics.
   - This step must precede read-support interpretation.

6. **Variant and coverage evidence**
   - Identify SNPs and indels relative to the current mitogenome.
   - Summarize depth, allele balance, support for insertions/deletions, and gene-level variant evidence.

7. **Read-support interpretation**
   - Interpret read mapping and variant evidence in the context of annotation problems.
   - Evaluate whether internal stops, frameshifts or suspicious CDS features are supported by reads or suggest sequence/annotation errors.

8. **Targeted consensus and local reconstruction**
   - Used for missing, partial or ambiguous genes.
   - Particularly relevant for genes such as a missing or poorly represented CDS.

9. **Candidate assembly**
   - Local assembly of candidate regions when read-level evidence and targeted consensus remain insufficient.

10. **Curated outputs and integrated report**
   - Produce final curated GenBank/FASTA outputs.
   - Generate final integrated curation reports documenting all evidence and actions.

This model prevents downstream read-support or candidate-assembly interpretation from being executed before the required read mapping and variant evidence layers exist.

---

## Current evidence-driven execution model

The MitoCurator workflow is organized as an evidence-driven curation pipeline. The intended order is:

1. **Input mitochondrial sequence**
   - FASTA sequence recovered by tools such as MitoHiFi, PipeASM or another assembly/recovery workflow; or
   - an already annotated GenBank file.

2. **Initial annotation**
   - If the input is FASTA, MitoFinder is used to generate an initial GenBank annotation.
   - If the input is GenBank, the existing annotation is used as the starting point.

3. **Initial refinement**
   - Standardize gene names and feature labels.
   - Compare the annotation against an expected mitochondrial gene set.
   - Add or normalize non-coding regions such as AT-rich/control-region features when supported.
   - Generate refinement candidates for missing genes or problematic CDS features.

4. **Diagnosis and annotation assessment**
   - Measure annotation completeness and feature-level quality.
   - Detect missing genes, problematic CDS features, internal stops, possible boundary issues, and missing tRNAs/rRNAs.
   - Generate a curator-facing annotation assessment report.
   - Outputs include:
     - `annotation_gene_inventory.tsv`
     - `annotation_review_targets.tsv`
     - `annotation_evidence_summary.tsv`
     - `annotation_assessment_report.md`
     - `annotation_assessment_report.html`

5. **Read mapping to the mitogenome**
   - Map available sequencing reads against the current mitochondrial reference.
   - Generate sorted/indexed BAM files, coverage summaries and mapping statistics.
   - This step must precede read-support interpretation.

6. **Variant and coverage evidence**
   - Identify SNPs and indels relative to the current mitogenome.
   - Summarize depth, allele balance, support for insertions/deletions, and gene-level variant evidence.

7. **Read-support interpretation**
   - Interpret read mapping and variant evidence in the context of annotation problems.
   - Evaluate whether internal stops, frameshifts or suspicious CDS features are supported by reads or suggest sequence/annotation errors.

8. **Targeted consensus and local reconstruction**
   - Used for missing, partial or ambiguous genes.
   - Particularly relevant for genes such as a missing or poorly represented CDS.

9. **Candidate assembly**
   - Local assembly of candidate regions when read-level evidence and targeted consensus remain insufficient.

10. **Curated outputs and integrated report**
   - Produce final curated GenBank/FASTA outputs.
   - Generate final integrated curation reports documenting all evidence and actions.

This model prevents downstream read-support or candidate-assembly interpretation from being executed before the required read mapping and variant evidence layers exist.

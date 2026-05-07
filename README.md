# MitoCurator

## Visão geral
MitoCurator é um pipeline para anotação, refinamento e diagnóstico de mitogenomas com foco curatorial (diagnóstico/manual review), evitando edição automática agressiva do GenBank.

Fluxo atual: **FASTA-first** com anotação inicial via MitoFinder, seguida de refinement, rotação e diagnóstico final.

## Workflow atual (`mitocurator run`)
1. `[1/5] Tool check`
2. `[2/5] MitoFinder annotation`
3. `[3/5] Annotation refinement`
4. `[4/5] Rotation`
5. `[5/5] Diagnosis`

## Principais diretórios de saída
- `00_logs/`
- `03_mitofinder/`
- `05_refinement/`
- `04_rotation/`
- `06_read_support/`
- `07_gene_qc/`

## Arquivos em `05_refinement/`
- `refined.gb`: GenBank refinado (diagnóstico-only; sem correção automática de genes/CDS).
- `expected_gene_set.tsv`: checagem do conjunto esperado.
- `added_features.tsv`: features adicionadas no refinement (ex.: AT-rich).
- `missing_gene_candidates.tsv`: candidatos para genes ausentes.
- `cds_refinement_candidates.tsv`: candidatos de refinamento para CDS problemáticas.
- `reference_similarity_candidates.tsv`: comparação por referência para candidatos de genes ausentes.
- `problematic_cds_reference_check.tsv`: comparação por referência para CDS com stop interno.
- `problematic_cds_stop_context.tsv`: contexto genômico dos códons stop internos.
- `problematic_cds_reference_alignment.tsv`: mapeamento de posições problemáticas vs referência.
- `missing_gene_candidate_proteins.faa`: proteínas traduzidas de candidatos de genes ausentes.
- `problematic_cds_proteins.faa`: proteínas de CDS problemáticas (mantendo `*`).
- `curation_recommendations.tsv`: tabela agregada de recomendações de curadoria por prioridade.
- `curation_recommendations.md`: relatório legível de recomendações e resumo final.

## Arquivos em `07_gene_qc/`
- `gene_qc.tsv`
- `problematic_features.tsv`
- `intergenic_regions.tsv`
- `diagnostic_summary.md`

## Configuração YAML (exemplos)
### a) Inseto/abelha
```yaml
project:
  genetic_code: 5

refinement:
  expected_gene_set:
    profile: metazoa_mito
    custom_file: null
  gene_name_profile: insect_mito
```

### b) Vertebrado/peixe/anfíbio
```yaml
project:
  genetic_code: 2

refinement:
  expected_gene_set:
    profile: vertebrate_mito
    custom_file: null
```

### c) Custom
```yaml
refinement:
  expected_gene_set:
    profile: custom
    custom_file: /path/to/expected_genes.tsv
```

### d) Referência GenBank
```yaml
reference:
  genbank: /path/to/reference.gb
```

## `expected_gene_set` profiles
- `metazoa_mito`
- `insect_mito`
- `vertebrate_mito`
- `minimal_mito`
- `custom`

> `expected_gene_set` **não** define o código genético.

## Código genético (ordem de prioridade)
1. `project.genetic_code`
2. `genetic_code`
3. `annotation.genetic_code`
4. `mitofinder.organism_code`
5. default `5`

## Comparação com referência
A comparação por referência:
- avalia candidatos de genes ausentes;
- avalia CDS problemáticas (ex.: stop interno);
- **não adiciona nem corrige genes automaticamente**;
- prioriza alvos para curadoria manual.

## Comando exemplo
```bash
python -m mitocurator.cli run --config config.teste.yaml
```

## Limitações atuais
- não corrige automaticamente o GenBank;
- não adiciona ND2 automaticamente;
- não corrige automaticamente CDS com stop interno;
- comparação com referência é diagnóstica;
- MitoFinder continua como dependência externa para anotação inicial.


## Arquivos em `06_read_support/`
- `problematic_stop_read_support.tsv`
- `problematic_stop_variants.tsv`
- `read_support_summary.md`

# MitoCurator

## Visão geral
MitoCurator é um pipeline para anotação, refinamento e diagnóstico de mitogenomas com foco curatorial (diagnóstico/manual review), evitando edição automática agressiva do GenBank.

Fluxo atual: **FASTA-first** com anotação inicial via MitoFinder, seguida de refinement, rotação e diagnóstico final.

## Workflow atual (`mitocurator run`)
1. `[1/5] Tool check`
2. `[2/5] MitoFinder annotation`
3. `[3/5] Annotation refinement`
4. `[4/5] Rotation`
5. `[5/6] Read support`
6. `[6/6] Diagnosis`

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


## Tecnologias de leitura suportadas em read_support
- PacBio HiFi
- PacBio CLR
- Oxford Nanopore (ONT)
- Illumina paired-end
- Illumina single-end

A configuração aceita formato legado (`use_hifi`/`use_illumina`) e formato novo multi-read-set (`read_support.read_sets`).
Não é necessário repetir caminhos dentro de `read_support`: use `source` para apontar para caminhos já definidos em `reads`.
O `refined.gb` é usado automaticamente pelo pipeline e convertido para `06_read_support/refined.fa` quando necessário.

### Exemplo recomendado (sem duplicar caminhos)
```yaml
reads:
  hifi:
    - /path/hifi.fastq.gz
  illumina:
    r1: /path/R1.fq.gz
    r2: /path/R2.fq.gz
  ont:
    - /path/ont.fastq.gz
  pacbio_clr:
    - /path/clr.fastq.gz
  illumina_se:
    - /path/illumina_se.fq.gz

read_support:
  enabled: true
  read_sets:
    - name: hifi
      type: pacbio_hifi
      source: reads.hifi
    - name: ont
      type: ont
      source: reads.ont
    - name: pacbio_clr
      type: pacbio_clr
      source: reads.pacbio_clr
    - name: illumina
      type: illumina_pe
      source: reads.illumina
    - name: illumina_se
      type: illumina_se
      source: reads.illumina_se
```


## Configuração simplificada (inputs no topo)
Defina no início do YAML: `project.output_base_dir`, `project.output_prefix`, `input.mitogenome`, `reference.genbank` e `reads.long`/`reads.short`.

- `reference.genbank` é reutilizado por MitoFinder/refinement/read_support quando necessário.
- `refined.gb` é gerado internamente e não precisa ser informado manualmente.
- Reads longas/curtas são declaradas uma única vez e reaproveitadas em `read_support.read_sets` via `source`.

## Novos arquivos em `06_read_support/`
- `readset_consensus_recommendations.tsv`: consenso por stop entre múltiplos read sets.
- `readset_consensus_recommendations.md`: relatório consolidado com conflitos e consenso.

Categorias de consenso incluem também suporte parcial entre read sets:
- `STOP_SUPPORTED_BY_SOME_READSETS`
- `CORRECTION_SUPPORTED_BY_SOME_READSETS`

## Cache de BAM na etapa read_support
Para acelerar testes iterativos, `read_support` aceita:

```yaml
read_support:
  reuse_existing_bam: true
  force_remap: false
```

- Se `reuse_existing_bam: true` e `{name}_to_refined.bam` + `.bai` já existirem, o pipeline reaproveita o BAM.
- Se `force_remap: true`, o mapeamento é refeito mesmo com BAM existente.

## Etapa `08_targeted_extraction/`
Nova etapa opcional para extrair reads direcionadas por alvo (diagnóstico/curadoria), usando **BAMs já produzidos em `06_read_support/`** (sem remapear reads e sem corrigir GenBank automaticamente).

Arquivos:
- `08_targeted_extraction/targets.bed`
- `08_targeted_extraction/targeted_read_extraction.tsv`
- `08_targeted_extraction/targeted_read_extraction.md`
- `08_targeted_extraction/reads/*.fastq.gz`

Para Illumina paired-end, a etapa exporta preferencialmente `*_R1.fastq.gz` e `*_R2.fastq.gz`; quando não há pareamento confiável suficiente, usa fallback `*.interleaved.fastq.gz` (registrado no TSV).

Tipos de alvo:
- candidatos de genes ausentes (`reference_similarity_candidates.tsv`, PARTIAL/STRONG),
- CDS problemáticas (`problematic_cds_stop_context.tsv`),
- candidatos de correção por consenso (`readset_consensus_recommendations.tsv`).

## Etapa `09_reconstruction_pools/`
Etapa opcional para montar **pools de reads** por alvo para remontagem/refinamento posterior (sem executar montadores automaticamente).

Entradas principais:
- `08_targeted_extraction/targeted_read_extraction.tsv`
- `06_read_support/{read_set}_to_refined.bam`

Pools gerados por alvo/read_set:
1. `target_only` (reads do alvo já extraídas),
2. `mitogenome_mapped` (todas as reads mapeadas no mitogenoma refinado),
3. `combined` (união deduplicada target_only + mitogenome_mapped).

Saídas:
- `09_reconstruction_pools/reconstruction_pools.tsv`
- `09_reconstruction_pools/reconstruction_pools.md`

> Nenhum montador é executado nesta etapa e nenhuma correção automática é aplicada ao GenBank.

## Etapa `10_targeted_consensus/`
Etapa opcional para construir consenso local por alvo/read_set (por padrão usando `pool_type: combined`) para apoiar decisão curatorial.

Saídas:
- `10_targeted_consensus/targeted_consensus.tsv`
- `10_targeted_consensus/targeted_consensus.md`
- `10_targeted_consensus/consensus_fasta/*.consensus.fasta`

Implementação atual (método `pileup`) usa maioria por posição com filtros de qualidade/profundidade e gera recomendações diagnósticas, sem alterar o GenBank.

### Rodando apenas targeted_consensus
```bash
python -m mitocurator.cli targeted-consensus --config config.yaml
```

### Teste rápido (ex.: apenas ND2 em hifi)
```yaml
targeted_consensus:
  target_filter: ND2
  read_set_filter: hifi
  max_targets: 2
```

### Presets minimap2 esperados por tecnologia
- HiFi (`pacbio_hifi`/`hifi`): `map-hifi`
- PacBio CLR (`pacbio_clr`/`clr`): `map-pb`
- ONT (`ont`/`nanopore`): `map-ont`
- Illumina PE/SE: `sr`

Parâmetros de desempenho:
- `jobs` (reservado para paralelização futura),
- `minimap2_threads`,
- `samtools_threads`.

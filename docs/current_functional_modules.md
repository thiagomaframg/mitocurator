# MitoCurator — módulos funcionais atuais (auditoria)

## Visão geral
Este documento descreve **apenas o que já está funcional no snapshot atual** da branch `codex/01-roadmap-and-functional-audit`.

## CLI e orquestração

### `mitocurator/cli.py`
Comandos disponíveis:
- `check-tools`: valida disponibilidade/configuração de ferramentas externas.
- `diagnose`: executa diagnóstico gene-a-gene e regiões intergênicas.
- `rotate`: rotaciona GenBank anotado para gene inicial definido pelo usuário.
- `run`: fluxo all-in-one inicial (check-tools → refinement opcional → rotate com tolerância a falha → diagnose).

Saídas principais esperadas em `run`:
- `00_logs/tool_check.tsv`
- `05_refinement/refined.gb`
- `05_refinement/*.tsv` de refinamento
- `07_gene_qc/gene_qc.tsv`
- `07_gene_qc/problematic_features.tsv`
- `07_gene_qc/intergenic_regions.tsv`
- `07_gene_qc/diagnostic_summary.md`

## Diagnóstico de features

### `mitocurator/gene_qc.py`
Responsabilidades principais:
- Ler registro (FASTA/GenBank) via camada de I/O.
- Resumir features não-`source` (tipo, coordenadas, strand, tamanho).
- Para `CDS`:
  - traduzir com código genético definido no config;
  - detectar stops internos e terminal;
  - marcar problema quando:
    - tamanho não é múltiplo de 3;
    - há stop interno;
    - ou ambos.
- Gerar regiões intergênicas lineares (entre features ordenadas), com tamanho e %AT.
- Escrever:
  - `gene_qc.tsv`
  - `problematic_features.tsv`
  - `intergenic_regions.tsv`
  - `diagnostic_summary.md`

## Refinamento de anotação (heurístico, assistivo)

### `mitocurator/refinement.py`
Responsabilidades observadas:
- Definir conjunto esperado de genes mitocondriais (CDS/rRNA/tRNA).
- Normalizar aliases de nomes de genes.
- Contar presença/duplicação/ausência em `expected_gene_set.tsv`.
- Inspecionar lacunas intergênicas e, se critérios forem atendidos, adicionar `misc_feature` de região AT-rich putativa.
- Buscar candidatos para CDS ausentes em lacunas intergênicas por varredura de ORFs (sem “resgate automático”).
- Buscar candidatos de ajuste para CDS com stops internos em janela local.

Importante:
- O módulo é de **suporte à decisão** e não substitui validação manual final.
- As tabelas de candidatos carregam `decision_hint`, servindo como priorização curatorial.

## Módulo de rotação

### `mitocurator/rotate.py`
- Rotaciona genoma circular anotado com base no gene inicial configurado.
- No fluxo `run`, falha de rotação não interrompe diagnóstico (pipeline segue adiante).

## Camadas de suporte
- `mitocurator/check_tools.py`: checagem operacional de ferramentas externas configuradas.
- `mitocurator/io.py`: leitura/escrita de registros e utilidades de nomenclatura.
- `mitocurator/utils.py`: utilitários de config e filesystem.

## Lacunas em relação ao plano expandido
Os módulos abaixo foram citados no objetivo macro, mas **não existem no estado atual do repositório**:
- `read_support.py`
- `reconstruction_pools.py`
- `targeted_extraction.py`
- `targeted_consensus.py`
- `candidate_assembly.py`

Conclusão: a pipeline atual é majoritariamente diagnóstica/refinamento leve; validações por suporte de leitura e reconstruções direcionadas ainda não estão implementadas neste snapshot.

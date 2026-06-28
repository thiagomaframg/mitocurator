# MitoCurator — Dev Brief: módulo `local_consensus`

Referência de implementação obrigatória para `mitocurator/local_consensus.py`.
Ler antes de qualquer commit no `feature/local-consensus-core`.

---

## Problema que resolve

O pipeline de refinement (`mitocurator/refinement.py`) diagnostica dois tipos de
problema em CDS:

1. **PROBLEM_INTERNAL_STOP** — CDS presente na anotação com ≥ 1 stop codon interno
   ao ser traduzida com o código genético do projeto.
2. **MISSING** — gene esperado ausente da anotação, ou presente com comprimento ou
   cobertura de aminoácidos insuficiente vs. a referência.

O módulo `local_consensus.py` resolve esses problemas usando evidência direta das
reads: extrai as reads que cobrem a região, constrói um alinhamento múltiplo com
MAFFT e deriva um consenso por maioria de coluna. Não aplica VCF global.

A motivação para MAFFT (e não pileup puro): pileup ancorado em coordenadas de
referência (`count_coverage`) resolve substituições, mas tem ponto cego para
frameshifts — um base a mais no assembly faz com que reads mapeiem com deleção
naquela posição, que `count_coverage` ignora, perpetuando o erro. MAFFT alinha
as reads como sequências completas e expõe colunas onde a maioria tem gap
(→ deletar base do assembly) ou onde a maioria tem base e o assembly tem lacuna
(→ inserir base). Foi a abordagem que funcionou na curadoria manual de
*Melipona capixaba*.

---

## Critério de acionamento (quando enviar um CDS para este módulo)

Um CDS é incluído na lista `problems` se satisfizer **qualquer** das condições:

| Condição | Como detectar |
|----------|---------------|
| Stop interno na tradução | `_translate_stop_metrics()` → `internal_stop_count ≥ 1` |
| Gene ausente da anotação | `summarize_expected_gene_set()` → `status == "MISSING"` |
| Comprimento < 80 % da ref | `len(CDS_nt) / ref_nt_len < completeness_min_ratio` (default 0.80) |
| Cobertura de aa < 70 % | identidade de alinhamento tblastn < `completeness_min_cov` (default 0.70) |

Os dois últimos critérios requerem proteína de referência em
`config.input.reference_gb`. Sem proteína de referência para o gene, registrar
`action: SKIPPED_NO_REFERENCE` no audit log e pular.

---

## Fluxo dentro do módulo (sem TSVs intermediários)

```
repair_cds_local_consensus(config, record, problems, reads_cfg, outdir, mode)
│
│  para cada CDS em problems:
│
├── [1] Extrair região ± flank_bp do mitogenoma → region.fa (FASTA temporário)
│         Incluir flancos para garantir contexto de mapeamento.
│
├── [2] minimap2 (preset por tecnologia) → samtools sort/index → BAM filtrado
│         Filtros: MQ ≥ min_mq, sem secondary (flag 256), sem supplementary (2048)
│
├── [3a] pysam.fetch() na região → recortar cada read por CIGAR/get_aligned_pairs()
│          para extrair apenas o segmento que se alinha às coordenadas da região
│          (incluindo flank_bp) → escrever segmentos recortados em reads_subset.fa.
│          Reads HiFi têm 10–20 kb; sem recorte o MAFFT alinharia sequências quase
│          inteiras contra uma region.fa pequena (lento e prejudicial ao alinhamento
│          da região de interesse). O recorte por aligned_pairs replica o que a
│          curadoria manual fez: extrair "a região ND1 de cada read", não a read
│          completa. Registrar no audit log: reads_trim_method="aligned_pairs".
│          Subsample aleatório até mafft_max_reads após o recorte.
│
├── [3b] mafft --auto --thread N --quiet <(cat region.fa reads_subset.fa)
│          → MSA em FASTA (region como primeira sequência)
│
├── [3c] Consenso coluna-a-coluna por maioria simples:
│          - Ignorar colunas onde fração de gaps > mafft_gap_col_threshold (0.50)
│            → essas colunas representam inserções no assembly (deletar base)
│          - Colunas com gap na maioria mas base em ≥ 1 read: inserir base majoritária
│          - Posições com profundidade < min_depth → 'N'
│          - Base chamada se freq ≥ maj_freq (default 0.70), senão 'N'
│
├── [4] Traduzir consenso com genetic_code do config
│         → contar stops internos
│         → calcular identidade aa vs. proteína de referência (PairwiseAligner)
│         → rejeitar candidato se stops internos > 0 ou N_fraction > max_n_fraction
│
└── [5] Audit log JSONL + saída conforme mode:
          diagnose → só audit log, sem arquivo de candidato
          suggest  → audit log + candidato em FASTA + anotação GenBank como
                     misc_feature com note="CANDIDATE_local_consensus"
          apply    → substitui coordenadas e sequência no GenBank in-place
                     (requer approved=true no dict do problema ou flag --apply)
```

---

## Assinatura pública

```python
def repair_cds_local_consensus(
    config: dict,
    record: SeqRecord,        # mitogenoma (modificado in-place se mode="apply")
    problems: list[dict],     # saída de find_cds_refinement_candidates() e
                              # find_missing_cds_candidates() do refinement
    outdir: Path,
    mode: str = "suggest",    # "diagnose" | "suggest" | "apply"
    audit_log: Path | None = None,
) -> list[dict]:              # entradas adicionadas ao audit log nesta chamada
```

`problems` usa os campos já produzidos pelo refinement:
`gene`, `seqid`, `start`, `end`, `strand`, `internal_stop_count`, `candidate_id`.
O módulo não re-executa o diagnóstico.

---

## Parâmetros de config (`config.yaml`)

```yaml
local_consensus:
  flank_bp: 200                  # flancos adicionados à região do gene
  min_mq: 20                     # mapping quality mínimo (samtools view -q)
  min_depth: 5                   # profundidade mínima para chamar base; abaixo → N
  maj_freq: 0.70                 # frequência mínima do alelo majoritário por coluna
  mafft_max_reads: 500           # máximo de reads passadas ao MAFFT (subsample)
  mafft_gap_col_threshold: 0.50  # fração de gaps para deletar coluna (frameshift)
  max_n_fraction: 0.05           # fração máxima de Ns no candidato para aceitar
  completeness_min_ratio: 0.80   # comprimento mínimo vs. referência
  completeness_min_cov: 0.70     # cobertura de aa mínima no alinhamento vs. ref
  threads: 8

reads:
  # Ao menos um dos grupos abaixo deve estar presente:
  hifi: /caminho/para/reads.fastq.gz      # caminho absoluto, ou lista de caminhos
  # illumina_r1: /caminho/R1.fastq.gz
  # illumina_r2: /caminho/R2.fastq.gz
  technology: hifi   # hifi | ont | clr | illumina
```

Preset minimap2 derivado de `technology`:

| technology | preset minimap2 |
|------------|-----------------|
| hifi       | `-x map-hifi`   |
| ont        | `-x map-ont`    |
| clr        | `-x map-pb`     |
| illumina   | `-ax sr`        |

---

## Estrutura de saída em `outdir`

```
06_local_consensus/
├── audit_log.jsonl                       ← append por execução
├── CYTB/
│   ├── region.fa                         ← região extraída do mitogenoma
│   ├── reads_subset.fa                   ← reads selecionadas (subsampled)
│   ├── mapped.bam / mapped.bam.bai
│   ├── mafft_aligned.fa                  ← MSA (region + reads)
│   ├── CYTB_candidate.fa                 ← consenso final (apenas em suggest/apply)
│   └── CYTB_candidate_protein.fa         ← tradução do candidato
├── ND2/
│   └── ...
└── summary.tsv                           ← uma linha por CDS, com resultado
```

O número do diretório (`06_`) deve ser configurável via `config.output.step_dirs`
para não assumir posição fixa no pipeline.

---

## Audit log JSONL

Arquivo: `{outdir}/06_local_consensus/audit_log.jsonl`
Modo: append (nunca sobrescrever — o log é cumulativo entre execuções).
Uma linha JSON por CDS processado:

```json
{
  "timestamp": "2026-06-28T12:34:56",
  "gene": "CYTB",
  "problem": "PROBLEM_INTERNAL_STOP",
  "evidence": {
    "type": "mafft_majority_consensus",
    "reads_technology": "hifi",
    "reads_source": "/caminho/Mcap_INMA.trimmed.fastq.gz",
    "region": "atg005843l_path_rc_rotated:11900-13200",
    "flank_bp": 200,
    "min_mq": 20,
    "reads_trim_method": "aligned_pairs",
    "reads_fetched": 312,
    "reads_used_in_mafft": 312,
    "mafft_gap_col_threshold": 0.50,
    "cols_deleted_as_insertion": 1,
    "cols_inserted_from_reads": 0,
    "consensus_n_fraction": 0.01,
    "stop_codons_before": 1,
    "stop_codons_after": 0,
    "ref_protein_identity_pct": 98.7,
    "ref_protein_coverage_pct": 100.0
  },
  "candidate": {
    "sequence_nt": "ATG...",
    "length_nt": 1143,
    "length_aa": 381,
    "internal_stops": 0,
    "terminal_stop": "yes",
    "fasta_path": "06_local_consensus/CYTB/CYTB_candidate.fa"
  },
  "action": "SUGGEST",
  "tools": {
    "minimap2": "2.28-r1209",
    "samtools": "1.21",
    "mafft": "7.526",
    "pysam": "0.22.1",
    "biopython": "1.84"
  },
  "commands": [
    "minimap2 -t 8 -x map-hifi region.fa reads.fastq.gz | samtools sort -o mapped.bam",
    "samtools view -b -q 20 -F 2308 mapped.bam -o mapped.filtered.bam",
    "samtools index mapped.filtered.bam",
    "mafft --auto --thread 8 --quiet input.fa > mafft_aligned.fa"
  ],
  "mitocurator_version": "0.1.0-dev"
}
```

Campos obrigatórios: `timestamp`, `gene`, `problem`, `evidence.type`,
`evidence.stop_codons_before`, `evidence.stop_codons_after`, `candidate`,
`action`, `tools`, `commands`.

Valores de `action`: `SUGGEST`, `APPLIED`, `REJECTED_STOPS_REMAIN`,
`REJECTED_HIGH_N`, `SKIPPED_NO_REFERENCE`, `SKIPPED_NO_READS`.

---

## Referências de implementação (estudar, não copiar)

- **Loop de pileup posição-a-posição** (`inspect/pr36`, `targeted_consensus.py`,
  ~L289–350): mostra como usar `pysam.AlignmentFile.count_coverage()` e o
  fallback com `aln.pileup()`. Útil como referência para o passo [2] (filtragem
  do BAM) e para calcular `reads_fetched`.

- **Evidência por codon em reads individuais** (`inspect/pr36`, `read_support.py`,
  ~L280–352): mostra como percorrer reads com `aln.fetch()` +
  `read.get_aligned_pairs()` e reconstruir codons por posição de referência.
  Útil para calcular `stop_codons_before` como validação independente.

Não usar: `reconstruction_pools.py`, `candidate_assembly.py` (Flye/SPAdes),
backends `samtools mpileup` de `targeted_consensus.py`.

---

## Critérios de regressão — caso *Melipona capixaba*

### Arquivos de teste

| Papel | Caminho absoluto |
|-------|-----------------|
| Entrada (pré-correção) | `/home/thiagomafra/projetos/asf/capixaba/montagem/pipeasm_run1/pipeasm/results/Assembly/Mitogenome/Solo_Asm/final_mitogenome.gb` |
| Reads HiFi (trimmed) | `/home/thiagomafra/projetos/asf/capixaba/montagem/pipeasm_run1/pipeasm/results/Trimming_QC/HiFi/Mcap_INMA.trimmed.fastq.gz` |
| Referência (*M. scutellaris*) | `/home/thiagomafra/projetos/asf/scutellaris/Mscutellaris_mtDNA.gb` |
| Alvo de regressão | `/home/thiagomafra/projetos/asf/capixaba/montagem/pipeasm_run1/pipeasm/results/Assembly/Mitogenome/Solo_Asm/polimento_hifi_mitogenoma/busca_nd2_reads_hifi/Melipona_capixaba_mitogenome_curated_final_circular.gb` |

### Saída mínima esperada

A partir do arquivo pré-correção, após `mitocurator run` com `local_consensus`
ativo (mode=apply):

| Gene | nt   | aa  | Stops internos |
|------|------|-----|----------------|
| ND1  | 882  | 294 | 0              |
| CYTB | 1143 | 381 | 0              |
| ND4  | 1305 | 435 | 0              |
| ND2  | ≈982 | —   | 0              |

Mitogenoma final: **19.526 bp**, circular, **13 CDS** sem stops internos,
**2 rRNA**, **22 tRNA** (anotação MiTFi).

### Problemas no arquivo de entrada que o módulo deve corrigir

| Gene | Tipo de problema | Stops internos (entrada) |
|------|-----------------|--------------------------|
| COX1 | PROBLEM_INTERNAL_STOP | 2 |
| ND6  | PROBLEM_INTERNAL_STOP | 5 |
| CYTB | PROBLEM_INTERNAL_STOP | 1 |
| ND5  | PROBLEM_INTERNAL_STOP | 4 |
| ATP8 | MISSING | — |
| ND2  | MISSING | — |

ATP8 e ND2 estavam presentes nas reads mas ausentes no consenso inicial — o
módulo deve recuperá-los mapeando contra a sequência de referência de
*M. scutellaris* e construindo consenso local.

### Critério de aprovação do teste

```bash
# após mitocurator run --mode apply:
python3 - <<'EOF'
import sys
from Bio import SeqIO

EXPECTED_13  = ["ATP6","ATP8","COX1","COX2","COX3","CYTB",
                "ND1","ND2","ND3","ND4","ND4L","ND5","ND6"]
EXACT_LEN    = {"ND1": 882,  "CYTB": 1143, "ND4": 1305}
EXACT_AA     = {"ND1": 294,  "CYTB": 381,  "ND4": 435}
GENETIC_CODE = 5  # ler de config na prática

rec = next(SeqIO.parse("output/final.gb", "genbank"))
cds = {}
for f in rec.features:
    if f.type == "CDS":
        name = f.qualifiers.get("gene", f.qualifiers.get("product", ["?"]))[0]
        cds[name] = f

failures = []

# 1. todos os 13 CDS presentes e sem stops internos
for gene in EXPECTED_13:
    if gene not in cds:
        failures.append(f"MISSING: {gene}")
        print(f"FAIL  {gene}: ausente")
        continue
    nt = cds[gene].extract(rec.seq)
    aa = str(nt.translate(table=GENETIC_CODE, to_stop=False))
    internal = aa[:-1].count("*")
    if internal:
        failures.append(f"STOPS: {gene} ({internal} stops internos)")
    print(f"{'OK  ' if not internal else 'FAIL'} {gene}: {len(nt)} nt / {len(aa)-1} aa / {internal} stops")

# 2. comprimento exato para genes de regressão
for gene in [g for g in EXACT_LEN if g in cds]:
    nt = cds[gene].extract(rec.seq)
    aa = str(nt.translate(table=GENETIC_CODE, to_stop=False))
    ok = len(nt) == EXACT_LEN[gene] and len(aa)-1 == EXACT_AA[gene]
    if not ok:
        failures.append(f"LEN: {gene} ({len(nt)} nt / {len(aa)-1} aa, esperado {EXACT_LEN[gene]} nt / {EXACT_AA[gene]} aa)")
    print(f"{'OK  ' if ok else 'FAIL'} {gene} comprimento: {len(nt)} nt / {len(aa)-1} aa")

# 3. comprimento do mitogenoma
# Ground truth específico deste caso de validação (conhecido pela curadoria manual prévia).
# O pipeline de produção NUNCA deve usar comprimento total esperado como critério de sucesso
# — espécies novas não têm essa resposta de antemão. Este check existe só para confirmar
# que o pipeline reproduz a resposta já conhecida deste caso.
ok_genome = len(rec.seq) == 19526
if not ok_genome:
    failures.append(f"GENOME_LEN: {len(rec.seq)} bp (esperado 19526)")
print(f"\nMitogenoma: {len(rec.seq)} bp {'OK' if ok_genome else 'FAIL (esperado 19526)'}")

print(f"\n{'REGRESSÃO OK' if not failures else 'REGRESSÃO FALHOU'}")
for f in failures:
    print(f"  - {f}")
sys.exit(len(failures))
EOF
```

---

## O que este módulo NÃO faz

- Não aplica correções globais de VCF.
- Não roda Flye, SPAdes ou qualquer montador de novo.
- Não altera o GenBank sem registrar no audit log.
- Não altera o GenBank em `mode="suggest"` (default).
- Não opera em `mode="apply"` sem flag explícito do usuário.
- Não infere código genético — usa sempre `config.project.genetic_code`.
- Não hardcoda posições, genes ou espécies específicas.

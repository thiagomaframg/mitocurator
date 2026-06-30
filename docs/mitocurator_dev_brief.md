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
  min_read_pident: 75.0          # piso de identidade por read no tblastn-vs-reads
                                 # (usado quando bimodalidade não é detectada)
                                 # Justificativa: 75% ~ divergência de família em proteínas
                                 # mitocondriais; reads abaixo disso vs. referência do mesmo
                                 # gênero são improváveis de ser mito real ou NUMT recente.
  min_bimodal_gap: 1.0           # gap mínimo (pp) para declarar bimodalidade na
                                 # distribuição de pident e cortar automaticamente no vale
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

Valores de `action`:

| Valor | Significado |
|-------|-------------|
| `APPLIED` | Candidato aceito e aplicado ao GenBank in-place |
| `SUGGEST` | Candidato aceito (mode=suggest), não aplicado |
| `REJECTED_STOPS_REMAIN` | Consenso final ainda contém stops internos |
| `REJECTED_HIGH_N` | Fração de Ns no candidato > `max_n_fraction` |
| `REJECTED_LOW_IDENTITY` | Identidade proteica vs. referência < `min_identity_pct` |
| `REJECTED_NO_BOUNDARY_FIT` | Nenhuma janela no boundary scan satisfaz os critérios |
| `REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY` | tblastn contra o assembly não confirma coordenadas |
| `SKIPPED_NO_REFERENCE` | Gene ausente na proteína de referência |
| `SKIPPED_NO_READS` | Nenhum arquivo de reads configurado |
| `SKIPPED_NO_READS_IN_BAM` | BAM filtrado contém 0 reads na região |
| `SKIPPED_NO_TBLASTN_HIT_IN_READS` | Nenhuma read passou o filtro de pident |
| `SKIPPED_NO_READS_PASS_FILTER` | Nenhuma read passou o filtro de stops por frame |
| `SKIPPED_NO_CONSENSUS` | MAFFT não produziu consenso válido |

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

| Gene | nt   | aa  | Stops internos | Status atual |
|------|------|-----|----------------|--------------|
| ND6  | 522  | 173 | 0              | APPLIED ✓   |
| CYTB | 1143 | 380 | 0              | APPLIED ✓   |
| COX1 | 1560 | 519 | 0              | APPLIED ✓ (anotação truncada — 702 nt → 1560 nt) |
| COX3 | 777  | 258 | 0              | APPLIED ✓ (anotação truncada — 477 nt → 777 nt) |
| ND1  | 882  | 293 | 0              | APPLIED ✓ (start codon canônico preferido — ATT/ATA ignorados) |
| ND4  | 1305 | 434 | 0              | APPLIED ✓ (filtro NUMT bimodal: 797 reads @88% descartadas) |
| ND5  | 1617 | 538 | 0              | APPLIED ✓ (anotação truncada) |
| ND2  | —    | —   | —              | REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY (M. scutellaris ND2 divergente demais) |

> **Nota ND4:** o aa esperado é 434 (*M. capixaba*), não 435 — a referência *M. scutellaris*
> tem 435 aa. O pipeline usa o comprimento da referência apenas para localizar a região e
> guiar o boundary scan; o comprimento final reflete a sequência real da espécie em curadoria.

Mitogenoma final: **19.526 bp**, circular, **13 CDS** sem stops internos,
**2 rRNA**, **22 tRNA** (anotação MiTFi).

> O comprimento 19.526 bp é específico do caso *M. capixaba* e não deve ser usado como
> critério de sucesso em espécies novas.

### Problemas no arquivo de entrada que o módulo deve corrigir

| Gene | Tipo de problema | Stops internos (entrada) | Resultado |
|------|-----------------|--------------------------|-----------|
| COX1 | INTERNAL_STOP | 2 | APPLIED |
| ND6  | INTERNAL_STOP | 5 | APPLIED |
| CYTB | INTERNAL_STOP | 1 | APPLIED |
| ND5  | INTERNAL_STOP | 4 | APPLIED |
| ND3  | INCOMPLETE_LENGTH | 0 | REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY (*M. scutellaris* ND3 divergente) |
| COX3 | INCOMPLETE_LENGTH | 0 | APPLIED |
| ND1  | INCOMPLETE_LENGTH | 0 | APPLIED |
| ND4  | INCOMPLETE_LENGTH | 0 | APPLIED (ver nota NUMT abaixo) |
| ATP8 | MISSING | — | não tentado (insect ATP8 ~168 nt, sobrepõe ATP6; ORF scan não encontra) |
| ND2  | MISSING | — | REJECTED_TBLASTN_NO_HIT_IN_ASSEMBLY (*M. scutellaris* ND2 divergente demais) |

#### Caso especial: ND4 e contaminação por NUMT

No dataset *M. capixaba*, a região do ND4 contém leituras de **dois grupos**:

- **Grupo NUMT** (~797 reads, pident médio 87.9% vs. *M. scutellaris*): cópia nuclear
  divergida; internamente consistente mas com ~9% de divergência em relação ao mito real.
- **Grupo mito** (~82 reads, pident ~98%): reads do mitogenoma real.

Sem filtro, a mistura produz 79 posições com frequência ~50/50 → Ns no consenso →
`REJECTED_HIGH_N`. A função `_find_pident_cutoff()` detecta a bimodalidade (gap de 1.05%
entre 93.9% e 94.9%) e corta automaticamente em 94.4%, descartando o cluster NUMT.
O resultado com apenas as 82 reads mitocondriais: 1305 nt / 434 aa / 0 stops / 0 Ns.

### Critério de aprovação do teste

```bash
# após mitocurator run --mode apply:
python3 - <<'EOF'
import sys
from Bio import SeqIO

EXPECTED_13  = ["ATP6","ATP8","COX1","COX2","COX3","CYTB",
                "ND1","ND2","ND3","ND4","ND4L","ND5","ND6"]
EXACT_LEN    = {"ND1": 882,  "CYTB": 1143, "ND4": 1305}
EXACT_AA     = {"ND1": 293,  "CYTB": 380,  "ND4": 434}  # valores M. capixaba, não M. scutellaris
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

---

## Detecção da região controladora (`refinement.py → add_at_rich_region`)

### Estratégia (três sinais em prioridade decrescente)

1. **Flancos de tRNA** (sinal primário, prior biológico). Parâmetro de config
   `refinement.control_region_flanks`: lista de pares `[tRNA_5prime, tRNA_3prime]`
   tentados em ordem. Se um par localiza o gap, o resultado é `confidence=high`
   independentemente do AT%.

   | Grupo taxonômico | Par canônico de flancos |
   |-----------------|------------------------|
   | Hymenoptera (default) | `["tRNA-Ile", "tRNA-Gln"]` |
   | Vertebrados (D-loop) | `["tRNA-Pro", "tRNA-Phe"]` |
   | Outros metazoários | consultar literatura; ajustar no config |

2. **Tandem repeat** (sinal secundário). Entre gaps candidatos com AT% ≥ threshold e
   comprimento ≥ min_len, preferir gaps com repetições em tandem detectáveis
   (~55 bp/unidade, ≥ 3 cópias, ≥ 75% identidade). `confidence=high` quando repeat
   detectado, `confidence=low` (revisar manualmente) quando apenas AT%.

3. **AT% isolado** (fallback final). `confidence=low`.

### Caso de validação — *Melipona capixaba* (HiFi, Hymenoptera)

- Entrada: `final_mitogenome.gb` (MitoFinder, pré-correção), 19506 bp.
- Anchor disponível: tRNA-Gln em 12104..12172. tRNA-Ile ausente (MitoFinder
  não anota tRNA-Ile em insetos; MiTFi necessário para anotação completa de tRNA).
- Gaps adjacentes a tRNA-Gln: 12173..14308 (2136 bp, 77.7% AT) e 12075..12103
  (29 bp, gap técnico entre rrnS e tRNA-Gln).
- Resultado: região 12173..14308 detectada como `confidence=high` pelo critério
  do maior gap adjacente.
- Falso positivo eliminado: gap 4006..6309 (antisense de COX1/ND4, gerado por
  anotação incorreta do MitoFinder que `local_consensus` corrige depois) não aparece
  como candidato porque o fallback AT+repeat não é executado quando o flanco localiza
  a região.

### Limitação conhecida — anchor único com gaps de tamanho similar

O critério "maior gap adjacente" é heurístico e foi validado apenas no caso
*M. capixaba*, onde a diferença de tamanho entre os gaps é nítida (2136 vs. 29 bp).
**Caso futuro de risco:** se um genoma tiver dois ou mais gaps de tamanho similar
(ex.: ±200 bp entre si) adjacentes ao único anchor disponível, a seleção pelo maior
se torna arbitrária e pode resultar em detecção incorreta sem qualquer aviso.

O comportamento correto para esse cenário seria:
- Classificar o resultado como `confidence=low` (ou rejeitar), e
- Registrar a ambiguidade no audit log para revisão manual.

**Não implementado.** Registrado aqui como item de revisão para quando um caso
ambíguo for encontrado — não especular sobre implementação antes disso.

---

## O que este módulo NÃO faz

- Não aplica correções globais de VCF.
- Não roda Flye, SPAdes ou qualquer montador de novo.
- Não altera o GenBank sem registrar no audit log.
- Não altera o GenBank em `mode="suggest"` (default).
- Não opera em `mode="apply"` sem flag explícito do usuário.
- Não infere código genético — usa sempre `config.project.genetic_code`.
- Não hardcoda posições, genes ou espécies específicas.

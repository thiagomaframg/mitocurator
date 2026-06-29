# Scripts de curadoria manual — *Melipona capixaba*

Referência dos scripts e comandos encontrados em:
```
/home/thiagomafra/projetos/asf/capixaba/montagem/pipeasm_run1/pipeasm/results/Assembly/Mitogenome/Solo_Asm/
```

Esses scripts documentam o fluxo de curadoria manual que o MitoCurator tenta automatizar.
Use este arquivo para cruzar a implementação com o procedimento original.

---

## 1. Polimento global do mitogenoma (antes da anotação)

### 1a. HiFi (2 rounds)

**Script:** `run_polish_mitogenome_hifi.sh`

Fluxo:
1. `minimap2 -ax map-hifi` — mapeia reads HiFi contra o mitogenoma
2. `samtools view -F 2308 -q 20` — filtra MAPQ ≥ 20, remove unmapped/secondary/supplementary
3. `bcftools mpileup -Q 13 -d 1000000` → `bcftools call --ploidy 1` — variant calling haploid
4. Filtra VCF: `QUAL ≥ 20`, `DP ≥ 5`
5. `bcftools consensus` — aplica variantes ao assembly
6. Repete (round 2 usa saída do round 1 como referência)

Cobertura atingida: ~10.965× (round 1), ~20.581× (round 2)
Saída: `polimento_hifi_mitogenoma/round2.hifi_polished.fasta`
Log: `run_polish_mitogenome_hifi.log`

### 1b. Illumina PCR-free (2 rounds)

**Script:** `run_polish_mitogenome_illumina.sh`

Fluxo:
1. `bwa-mem2` (ou `bwa mem`) — mapeia R1/R2 paired-end
2. `samtools view -q 20` — filtra MAPQ ≥ 20
3. `bcftools mpileup` → `bcftools call --ploidy 1`
4. Filtra VCF: `QUAL ≥ 30`, `DP ≥ 10`
5. `bcftools consensus` (2 rounds)

Saída: `polimento_illumina_mitogenoma/round2.polished.fasta`
Log: `run_polish_mitogenome_illumina.log`

---

## 2. Filtragem de reads de alta confiança via tblastn de proteína vs reads

Localização: `polimento_hifi_mitogenoma/busca_nd2_reads_hifi/`

Fluxo completo (aplicado a CYTB e ND1; mesma lógica para ND2):

### Etapa 1 — tblastn da proteína de referência contra reads candidatas

Referência extraída de: `polimento_hifi_mitogenoma/Mscutellaris_mtDNA.gb`
Reads-alvo: lista de read IDs mitocondriais pré-filtrada (`mito_candidate_read_ids.txt`, 154.036 reads; `.very_strict.txt`, 5.044 reads)

Saídas intermediárias:
- `cytb_tblastn_vs_reads.tsv` — 1.051 hits brutos
- `cytb_highconf_hits.tsv` — 61 hits de alta confiança
- `nd1_tblastn_vs_reads.tsv` — 939 hits brutos
- `nd1_highconf_hits.tsv` — 21 hits de alta confiança
- `nd1_candidate_read_ids.len800_mapq30.txt` — 278 reads
- `nd2_candidate_read_ids.len500_mapq30.txt` — 88 reads

### Etapa 2 — Extração da região alinhada de cada read

**Script CYTB:** `extract_cytb_regions_from_reads.py`
**Script ND1:** `extract_nd1_regions_from_reads.py`

Lógica: parseia o TSV do tblastn (read_id, coordenadas sstart/send, strand), extrai a subsequência correspondente de cada read candidata (orientação codificante).

Saídas: `cytb_highconf_regions.fasta`, `nd1_highconf_regions.fasta`

### Etapa 3 — Tradução e filtragem por stops internos

**Script CYTB:** `translate_cytb_regions.py`
- Traduz com **código genético 5** (mitocondrial de invertebrados)
- Busca nos 3 frames de leitura
- Mantém reads com **≤ 1 stop interno** (permite 1 terminal)
- Colunas de saída: id, nt_len, frame, aa_len, n_internal, internal_stops, first10aa, last10aa

**Script ND1:** `translate_nd1_regions.py`
- Mesmo fluxo, filtra para **0 stops internos** (critério mais estrito)

### Etapa 4 — MAFFT + consenso por maioria

**Script:** `consensus_from_cytb_alignment.py` / `consensus_from_alignment.py`
- Entrada: MSA (MAFFT das regiões filtradas)
- Consenso coluna a coluna, base mais frequente (exclui gaps)
- Saída: única sequência consenso FASTA

---

## 3. Refinamento de fronteiras CDS por eliminação de stops

Localização: `.../busca_nd2_reads_hifi/polished_contig3_test/.../`
Scripts auxiliares: `scan_orfs_table5.py`, `translate_3frames_table5.py`, `locate_stop_codons.py`, `check_internal_stops_table5.py`

### ND1 — `refine_nd1_by_stops.py`

Coordenadas base (de tblastn): `start=17280, end=18163, strand=+`
Lógica:
- Testa todos os shifts ±60 bp em cada extremidade
- Para cada combinação, testa os 3 frames de leitura
- Aceita se: **0 stops internos** E comprimento proteico **280–320 aa**

### CYTB — `refine_cytb_by_stops.py`

Coordenadas base (de tblastn): `18232..19365, reverse strand`
Lógica:
- Shifts ±30 bp em cada extremidade
- Extrai RC para genes na fita negativa
- Aceita se: **≤ 1 stop interno** E comprimento proteico **360–390 aa**

---

## 4. Análise de indels por posição (auxiliar)

**Script:** `polimento_illumina_mitogenoma/summarize_indels_from_bam.py`

Parseia BAM com pysam, reporta por posição: depth, insertions, deletions, ins_freq, del_freq.
Usado para identificar regiões problemáticas antes do polimento.

---

## 5. Extração de genes de GenBank (auxiliares)

Todos em `polimento_hifi_mitogenoma/`:

- `extract_proteins_from_gb.py` — extrai traduções CDS → FASTA de proteínas
- `extract_gene_nt_from_gb.py` — extrai sequência nucleotídica CDS (match exato de nome)
- `extract_gene_from_gb.py` — extrai sequência nucleotídica CDS (match por substring)

---

## 6. Fluxo completo resumido

```
Assembly bruto
    ↓
[Polimento HiFi — 2 rounds]       run_polish_mitogenome_hifi.sh
    ↓
[Polimento Illumina — 2 rounds]   run_polish_mitogenome_illumina.sh
    ↓
[Anotação — MitoFinder + MiTFi]   run_mitofinder_hifi_mscutellaris_ref_local_arwen.sh
    ↓
[Verificação de stops]             check_internal_stops_table5.py
    ↓ (para cada gene com stop)
[tblastn proteína vs reads]        → cytb/nd1_tblastn_vs_reads.tsv
[Extração região nas reads]        extract_cytb/nd1_regions_from_reads.py
[Filtragem por 0 stops internos]   translate_cytb/nd1_regions.py
[MAFFT + consenso por maioria]     consensus_from_cytb_alignment.py
[Refinamento de fronteiras]        refine_cytb/nd1_by_stops.py
    ↓
[Verificação circularização]       circularization_check_tm/check_mito_circularization.sh
    ↓
Mitogenoma curado final:
    Melipona_capixaba_mitogenome_curated_final_circular.gb
```

---

## Notas de implementação para o MitoCurator

- **Polimento global (bcftools)** → etapa anterior ao `repair_cds_local_consensus`; já rodado manualmente no caso M. capixaba; pipeline deve receber assembly já polido como entrada por padrão.
- **tblastn vs reads** (não vs assembly) → estratégia diferente do MitoCurator atual, que usa tblastn vs assembly para coordenadas e depois reads vs CDS de referência. A abordagem manual vai mais longe: usa as próprias reads como substrato do consenso, descartando o assembly inteiramente para aquele gene.
- **Código genético 5** em todos os scripts de tradução — consistente com `project.genetic_code` do MitoCurator.
- **CYTB**: o script manual aceitou ≤1 stop (critério relaxado); o MitoCurator atual exige 0 stops (mais rigoroso, usa rescue por referência para atingir isso).

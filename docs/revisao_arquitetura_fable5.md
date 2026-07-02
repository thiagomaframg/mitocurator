# Revisão crítica de arquitetura — MitoCurator

> Revisão solicitada antes de expandir o pipeline para novos táxons (peixes Neotropicais,
> anfíbios Microhylidae/Hylidae) e finalizar o manuscrito (Molecular Ecology Resources,
> Resource Article). Foco: encontrar falhas reais de arquitetura, não validar.
>
> Método: cruzamento da descrição do projeto com o código real (`mitocurator/*.py`),
> os scripts da curadoria manual de *Melipona capixaba* (`.../Solo_Asm/` e subdiretórios)
> e as issues abertas (#52, #53, #55). O audit log real (`output_test/06_local_consensus/`)
> foi inspecionado para ancorar o Problema A em números.
>
> **Ressalva de honestidade:** nenhum PDF de curadoria manual foi localizado no servidor
> (busca em `mitocurator/` e `asf/capixaba/` retornou vazio). O método-ouro foi cruzado
> contra os scripts `.sh`/`.py` da curadoria, não contra o PDF.

---

## Resumo por severidade

| # | Severidade | Achado |
|---|-----------|--------|
| 1 | **CRÍTICO** | Problema A não é "comprimento correto do dataset": é um bug de definição de borda que **duplica** sequência codificante e produz molécula biologicamente errada, invisível ao QC. Ameaça a tese central do manuscrito. |
| 2 | **ALTO** | `gene_qc` não tem QC de coordenadas/sobreposição/duplicação — só traduz CDS anotado. O ponto cego do MitoFinder (Q6) **e** o rastro de duplicação do achado #1 passam despercebidos. |
| 3 | **ALTO** | `polish.py` (bcftools ploidy-1, sem filtro NUMT) pode **fixar o alelo NUMT** numa região onde NUMT domina em profundidade, corrompendo o assembly *antes* do `local_consensus` agir. |
| 4 | **ALTO** | Recuperação de MISSING por reads falha silenciosamente para referência distante (Q2): não distingue "gene ausente" de "referência divergente demais". |
| 5 | **ALTO** | tRNA por contagem, não cobertura (Q3/#55) — e o desenho ARWEN-primário assume 22 universal e ARWEN majoritariamente certo. |
| 6 | **MÉDIO** | Tratamento de origem/circularidade linear em vários módulos — provável causa dos 5 tRNAs "ausentes" e risco em gene na junção. |
| 7 | **MÉDIO** | Reprodutibilidade do audit log insuficiente para MolEcolRes (Q5): sem hash/accession da referência, sem md5 do input, sem config, versão hardcoded. |
| 8 | **MÉDIO** | Não-determinismo (MAFFT `--auto` multi-thread, desempate de `Counter`, ordem de reads): o `output_test` diverge do run do PR #51. Mina a alegação de "reprodutível". |
| 9 | **MENOR** | Divergências código↔dev_brief (subsample de reads não implementado; trim por `aligned_pairs` não implementado). |

---

## 1. CRÍTICO — Problema A: causa raiz é bug de borda, não comprimento legítimo (Q4)

**Cadeia causal, confirmada no código:**

1. `refinement.py:531-543` (`find_cds_refinement_candidates`, ramo de completude) emite genes
   truncados como `INCOMPLETE_LENGTH` com `candidate_start/end` = **coordenadas da anotação
   truncada do MitoFinder** (`delta_start=0, delta_end=0`) e o comentário
   `"tblastn step in local_consensus will refine coordinates"`.
2. Em `local_consensus.py:126`, o tipo é decidido só por stops:
   `prob_type = "PROBLEM_INTERNAL_STOP" if stops_bef > 0 else "MISSING"`. Como esses genes têm
   0 stops, viram **`MISSING`**.
3. O ramo `MISSING` (`local_consensus.py:147-152`) **pula o tblastn de re-ancoragem no assembly**
   e faz `coords_confirmed = True` incondicionalmente. Logo `start0/end0` permanecem os da
   anotação truncada.
4. `_apply_to_record` (linha 822) executa `record.seq[:start0] + consenso_full + record.seq[end0:]`.
   Excisa o span **curto** e insere o consenso de **comprimento-referência**.

O comentário do refinement é uma **promessa falsa**: o tblastn prometido nunca roda para esses genes.

**Prova aritmética:** a issue #52 decompõe o delta em ND4 (+582), ND1 (+234), COX3 (+300),
ND3 (+219), ND2 (+150). O dev_brief diz literalmente que COX3 é "anotação truncada — 477 nt →
777 nt". `777 − 477 = 300` = exatamente o delta de COX3. Não é coincidência: o delta de cada gene
é precisamente `comprimento_real − comprimento_truncado`.

**Por que é biologicamente errado, não só "21 kb":** o gene completo (777 nt de COX3) quase
certamente **já está fisicamente no assembly** — o MitoFinder só anotou o boundary curto (parou
num ATG interno ou num stop prematuro por erro de assembly). Ao excisar 477 e colar 777, o
pipeline **duplica os ~300 nt de COX3 que já existiam a jusante** do `end0`. O resultado tem:
(a) +300 bp de comprimento, e (b) um trecho de COX3 duplicado, não anotado, na região intergênica
seguinte. O `gene_qc` valida só o CDS anotado (agora perfeito, 777 nt/0 stops) e **não vê** a cauda
duplicada. Isso é falso-verde.

Contraste com o método manual (`replace_cytb_in_ND1fixed.py`): o curador fixou
`replace_start=18237, replace_end=19370` — **as coordenadas reais e completas do gene no contig**,
determinadas por inspeção. Consenso de 1134 nt substituiu região de 1134 nt → delta ≈ 0. O humano
garantiu que o span excisado batia com a extensão real do gene. O pipeline nunca faz essa
verificação para MISSING/INCOMPLETE.

**Por que o caminho INTERNAL_STOP não infla:** ali o tblastn no assembly re-ancora `start0/end0`
para o gene inteiro (`_best_tblastn_coords`), então `old_span ≈ full` e `delta ≈ 0` — batendo com
COX1/CYTB/ND5/ND6 no audit log. **A correção é unificar: rodar o tblastn de re-ancoragem no
assembly para TODOS os tipos de problema.** Só quando o tblastn realmente não achar nada (ausência
verdadeira, como o ND2 real) é que se insere — e aí a junção precisa ser explicitamente validada e
sinalizada, não colada às cegas nas coordenadas de uma ORF de scan.

**Agravantes secundários no mesmo módulo:**
- `_extract_hit_regions` (linhas 601-602) **estende** cada read por `(qstart-1)*3` e
  `(ref_prot_len-qend)*3` para atingir o comprimento da referência. O manual
  (`extract_cytb_best_regions_from_reads.py`) extraía **só o bloco alinhado** (`sstart–send`).
  A extensão do pipeline puxa flanco não-codificante para dentro do consenso e permite ultrapassar
  a borda real do gene.
- `_scan_boundary_shifts` aceita candidatos dentro de ±15% do comprimento da referência e
  **prefere** início ATG/TTG, mas **não exige terminar em stop**. Um candidato pode não ser um CDS
  completo (sem stop terminal), e o comprimento fica atado à referência (espécie divergente!) e não
  à biologia real da espécie curada.

**Impacto no manuscrito:** a alegação central é "reproduz a curadoria manual, evidence-driven".
Hoje o pipeline **não reproduz** o padrão-ouro (21 kb com duplicações vs. 19.526 bp). Isto precisa
ser resolvido antes de submeter, senão o próprio caso de validação refuta a ferramenta.

---

## 2. ALTO — `gene_qc` não valida coordenadas nem estrutura (Q6, Q8)

`gene_qc.summarize_features` só faz: tradução de CDS (stops, múltiplo de 3) e comprimento de tRNA
(50-200 nt). **Não há**:
- checagem de **sobreposição/duplicação** de features (não pegaria a cauda duplicada do achado #1);
- verificação de que os **13 CDS esperados** estão todos presentes e únicos (isso está no
  refinement, mas não é reafirmado no QC final);
- QC de **coordenadas independente do MitoFinder**.

Sobre o ponto cego do MitoFinder (Q6): sim, é inaceitável como está para um "Resource Article".
Se o MitoFinder desloca um gene em +3 nt sem gerar stop interno, o gene passa por todo o pipeline
intacto e sai errado. O pipeline **já tem** a ferramenta para fechar esse buraco — o tblastn da
proteína de referência contra o assembly (`_best_tblastn_coords`) — mas só o usa quando há stop
interno. **Recomendado:** rodar um QC de coordenadas por tblastn/proteína para **todos** os CDS no
`gene_qc` final, comparando start/end anotados com o hit da referência e flagrando desvios (e
sobreposições e duplicações intergênicas de conteúdo codificante). Isso mata os achados #1 e #6 de
uma vez.

---

## 3. ALTO — Fronteira polish × local_consensus e o risco NUMT no polish (Q1)

A fronteira conceitual está **no lugar certo**: polish = substituição/indel global pré-anotação
(bcftools ploidy-1); local_consensus = reparo estrutural gene-a-gene pós-anotação (MAFFT). O
`polish.py` espelha fielmente `run_polish_mitogenome_hifi.sh` (mesmos filtros, 2 rounds, mesma
cadeia mpileup→call→norm→filter→consensus) — boa aderência.

**Mas há uma lacuna perigosa, não um overlap.** O `polish` chama variante por **maioria haploide**
(`bcftools call --ploidy 1`) e **não tem filtro NUMT**. Ele só *registra* `ambiguous_sites` na
faixa 30-70% de AF — mas não os corrige nem os usa. No caso ND4, o grupo NUMT tem ~797 reads vs.
~82 mito (~90% NUMT). Numa razão dessas, o ploidy-1 **chama o alelo NUMT como consenso** e o AF
(~90%) nem cai na faixa 30-70%, então não é sequer sinalizado. Ou seja: **o polish pode ativamente
converter a região do ND4 para a sequência NUMT antes do local_consensus rodar**, e o filtro
bimodal do local_consensus opera sobre *reads*, não desfaz o dano já gravado no assembly.

No fluxo manual isso não mordeu porque o polish rodou sobre o **contig mito já extraído**
(`final_mitogenome.fasta`), onde a maioria das reads NUMT não mapeia bem. Se no pipeline
`polish.assembly_fasta` for o assembly bruto completo, o risco é real. **Recomendado:**
(a) documentar que polish deve rodar sobre o contig mito, não o assembly bruto; e/ou (b) trazer o
mesmo racional bimodal do local_consensus para o polish (ou ao menos flagrar regiões com desvio
sistemático de profundidade/identidade), tratando NUMT como ameaça de primeira classe **também** no
polish — hoje ela só é tratada no local_consensus.

---

## 4. ALTO — Recuperação de MISSING não distingue "ausente" de "referência divergente" (Q2)

O caminho MISSING depende de `tblastn` (proteína da referência) achar hits nas reads, com
`evalue=1e-5` fixo e piso `min_read_pident=75%`. Para anfíbios com referência de **família**
distante, dois modos de falha:
- **Falso negativo silencioso:** nenhuma read passa 75% → `REJECTED_NO_READ_EVIDENCE`, que o log
  descreve como "gene provavelmente ausente ou divergente demais". Mas o gene **existe** — foi só a
  referência que não alcançou. O usuário lê "ausente" e conclui perda gênica biológica falsa.
- **Coordenada errada aceita:** um hit fraco e parcial vira âncora e o candidato é montado sobre
  bordas ruins.

O problema de fundo: **não há um sinal explícito de "referência divergente demais"** separado de
"gene ausente". `min_identity_pct=70` também rejeita um gene real que legitimamente diverge >30% da
referência de família. **Recomendado:** tornar o limiar explícito e reportado — p.ex. calcular a
identidade proteica esperada referência↔alvo (a partir de genes que *deram* certo, como COX1) e, se
a divergência global for alta, (a) baixar automaticamente os pisos, (b) preferir uma referência
congênere/co-familiar quando existir, e (c) emitir um status distinto tipo
`INCONCLUSIVE_REFERENCE_TOO_DIVERGENT` em vez de `REJECTED_NO_READ_EVIDENCE`. Para o manuscrito com
peixes e anfíbios, esse é o ponto que mais provavelmente gera conclusão biológica incorreta.

---

## 5. ALTO — tRNA por contagem é sintoma de desenho frágil (Q3, #55)

A issue #55 está certa: trigger e completude por **contagem bruta** deixam passar "22 = 18 reais +
4 duplicatas". Corrigir para **cobertura por (aminoácido+anticódon)** é necessário. Mas o problema é
mais profundo:

- **`expected_trna_count=22` como constante universal** não generaliza. Metazoários têm
  perdas/ganhos de tRNA; alguns têm 23-24; muitos vertebrados têm 22 mas com identidades distintas.
  Fixar 22 e o set `EXPECTED_GENE_SET` (hardcoded em `refinement.py:17-21`, com `Leu2/Ser2`) é um
  **hardcode de conteúdo gênico** que colide com o princípio "nada de valores de espécie no código".
  Deveria vir do config por táxon.
- **ARWEN-primário + MiTFi-só-complementa** (`mitfi.py:341-345`): se o ARWEN **anota um aminoácido
  com identidade errada**, o MiTFi nunca corrige, porque ele só adiciona AAs *ausentes* do conjunto
  ARWEN. Para táxons onde o ARWEN tem alto FN **ou** troca identidades, essa hierarquia propaga o
  erro. O CLAUDE.md diz "preferir MiTFi"; o código faz o oposto (ARWEN manda). Vale reconsiderar
  inverter: MiTFi (covariância, mais específico) como primário.
- **Os 5 ausentes (Ala/Cys/Glu/Lys/Met)** provavelmente **não** são biológicos. Ver achado #6: o
  MiTFi recebe a molécula como **FASTA linear** (`SeqIO.write(record, ..., 'fasta')`, linha 320) —
  qualquer tRNA na junção circular é cortado e não detectado. Antes de concluir "faltam 5 tRNAs",
  teste com a sequência duplicada/rotacionada (o próprio projeto tem `duplicate_circular_fasta.py`)
  e com `-both`/modelos alternativos, como a própria #55 sugere. A junção é candidata forte a causa
  de pelo menos parte.

---

## 6. MÉDIO — Circularidade tratada de forma linear em vários pontos (Q8)

Ordem de execução: `mitfi → refinement → local_consensus → rotate`. Todas as operações de gene
acontecem **antes** da rotação, sobre a origem arbitrária do assembly. Pontos frágeis:
- `local_consensus._extract_region` e `_apply_to_record` são **puramente lineares** — um gene com
  `start0 > end0` (cruzando a origem) ou próximo dela seria mal-extraído e o splice corromperia a
  molécula. Gene na junção não é suportado.
- `mitfi.py`: FASTA linear → tRNAs na origem perdidos (provável causa parcial dos 5 ausentes).
- `_best_tblastn_coords` não representa hit que enrola na origem.
- `find_missing_cds_candidates` (refinement) usa `(gap["start0"]+i) % len` e pode gerar
  `start > end` silenciosamente para gaps que enrolam.

`refinement._intergenic_gaps` e `final_molecule_preparation.detect_at_rich_region` **tratam** wrap
corretamente — então o tratamento é inconsistente entre módulos. **Recomendado:** decidir
explicitamente a política (p.ex. rotacionar para um marco estável **antes** de refinement/
local_consensus, garantindo que nenhum gene-alvo fique na junção) e documentar que gene-na-origem
não é suportado até haver teste.

---

## 7. MÉDIO — Reprodutibilidade do audit log para MolEcolRes (Q5)

Registra versões de ferramentas e comandos — bom, mas **insuficiente** para o padrão de
reprodutibilidade da Wiley/MolEcolRes. Faltam:
- **Accession + md5/hash da referência** usada (hoje só decisões; a referência é o determinante
  silencioso de todas as coordenadas e limiares).
- **md5 do assembly de entrada** e das **reads**.
- **Snapshot do config** (todos os limiares) no log.
- **Versões de TODAS as ferramentas** por etapa: falta bcftools/bwa no log do local_consensus;
  falta java/**MiTFi**/**infernal** no log do mitfi; falta tblastn/makeblastdb.
- `mitocurator_version` está **hardcoded `"0.1.0-dev"`** (`local_consensus.py:973`,
  `polish.py:191`) — deveria ser o **git hash**.
- **Seed aleatória**: ver achado #8.

Sem accession+hash da referência, um revisor não consegue re-executar de forma determinística — e a
resposta muda com a referência.

---

## 8. MÉDIO — Não-determinismo (Q8, e mina a alegação "reprodutível")

Observação empírica: o audit log de `output_test` mostra deltas pequenos (~+130 bp), enquanto o run
do PR #51 deu +1.494 bp. Parte é config (INCOMPLETE vs MISSING via ORF-scan), mas há fontes
genuínas de não-determinismo:
- `mafft --auto --thread N` pode variar a árvore/alinhamento entre execuções multi-thread.
- `Counter(...).most_common(1)` (`local_consensus.py:690`) desempata por ordem de inserção →
  depende da ordem das reads.
- `_extract_reads_from_bam` despeja reads na ordem do BAM; sem seed.

Para um "Resource Article" que vende reprodutibilidade, isso precisa de: MAFFT com seed/`--thread 1`
no modo de validação, desempate determinístico explícito, e ordenação estável das reads. Idealmente
um teste de regressão que rode 2× e exija saída idêntica.

---

## 9. MENOR — Divergências código ↔ dev_brief e outros

- O dev_brief (§3a) descreve **recorte por `aligned_pairs`** e **subsample até
  `mafft_max_reads=500`** em `reads_subset.fa`. No código, `_extract_reads_from_bam` despeja
  **reads inteiras** (10-20 kb) sem subsample; o recorte real vem depois via tblastn +
  `_extract_hit_regions`. Funciona, mas: (a) sem subsample, cobertura alta → MAFFT com milhares de
  sequências pode estourar tempo/memória; (b) a doc não descreve o que o código faz. Alinhar código
  e brief.
- `ATP8` nunca é tentado (overlap com ATP6, ORF-scan não acha) — genes pequenos/sobrepostos são um
  buraco assumido; ok documentar, mas é limitação real para generalização.
- `_apply_to_record` usa `int(feat.location.start/end)` — para features com `CompoundLocation`
  (junção), isso colapsa para min/max e desloca errado.
- `check_tools`: confirmar cobertura de **java/MiTFi/infernal/bwa/bcftools** como dependências de
  primeira classe, já que o `run` depende delas.

---

## Prioridade recomendada antes de expandir táxons e submeter

1. **Corrigir o achado #1** (unificar a re-ancoragem por tblastn no assembly para MISSING/INCOMPLETE;
   só inserir sem âncora em ausência verdadeira, com junção validada). Pré-requisito para reproduzir
   19.526 bp e sustentar o manuscrito. *(recomendado)*
2. **Adicionar QC de coordenadas/sobreposição independente do MitoFinder** no `gene_qc` (fecha #1 e
   #6). *(recomendado)*
3. **Tornar explícito o limite de divergência da referência** (#4) — crítico para peixes/anfíbios.
4. **Testar a hipótese da junção circular** nos 5 tRNAs (#5/#6) antes de concluir perda biológica.
5. **Endurecer reprodutibilidade** (#7/#8): hash da referência, git hash, determinismo.
6. Reavaliar **NUMT no polish** (#3) e a hierarquia **ARWEN/MiTFi** (#5).

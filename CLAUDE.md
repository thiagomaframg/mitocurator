# CLAUDE.md — MitoCurator

Instruções de projeto para o Claude Code, lidas automaticamente no início de cada sessão.
Complementam (não substituem) as instruções globais do usuário.

## O que é

MitoCurator: pipeline de curadoria de genomas mitocondriais de metazoários montados por
long reads, short reads ou dados híbridos. Diferencial sobre MITOS2/MitoFinder/MitoZ:
não apenas anota — **repara** a sequência com base em evidência das reads e **registra
cada decisão** (curadoria "evidence-driven"). Origem: curadoria manual do mitogenoma de
*Melipona capixaba*.

## Decisões de arquitetura (fixas — não reabrir sem consultar o usuário)

1. **Operação por etapa em três níveis:** `diagnose` (só reporta) / `suggest` (gera
   candidato + evidência; default) / `apply` (aplica candidato já aprovado). Nunca alterar
   sequência ou anotação sem registrar evidência e, por padrão, sem confirmação.
2. **Audit log obrigatório (JSONL):** toda decisão de curadoria registra gene, problema,
   evidência (reads/alinhamento/tradução), candidato proposto, ação tomada, comandos e
   versões de ferramentas. É a espinha dorsal do argumento do manuscrito.
3. **Generalidade para metazoários.** Código genético sempre de `project.genetic_code`
   (5 = invertebrado, 2 = vertebrado). **Proibido hardcode** específico do caso Melipona
   (ex.: posições fixas, stop codons fixos — derivar de `Bio.Data.CodonTable`).
4. Antes de mexer no núcleo (mapeamento → consenso-local), **ler
   `docs/mitocurator_dev_brief.md`** — especificação detalhada e critérios de regressão.

## Regras científicas

- **Nunca aplicar VCF global cegamente** para corrigir CDS — o polimento global falhou em
  ND1/CYTB no caso real. O reparo é local, por reads de alta confiança (tblastn da proteína
  de referência → extração → MAFFT → consenso por maioria → validar 0 stops internos).
- Gene esperado **ausente** na anotação → buscar evidência nas reads antes de concluir
  ausência (caso ND2: existia nos dados, não no consenso inicial).
- **Não inventar dados:** accession numbers, nomes de genes/espécies, resultados ou
  referências. Sem certeza → dizer "não sei".
- **Nomenclatura:** nomes científicos em itálico na documentação; normalizar nomes de gene
  (COI/COX1, cob/CYTB, ND/NAD…) pela função existente. Em insetos usar
  "control region"/"A+T-rich region", **nunca "D-loop"**.

## Convenções de código

- **Simplicidade primeiro:** o mínimo de código que resolve. Nada especulativo; sem
  abstração para uso único; sem parâmetros/flexibilidade não pedidos; sem tratamento de
  erro para cenários que não ocorrem no fluxo real.
- **Mudanças cirúrgicas:** tocar só no necessário. Não "melhorar" código adjacente, não
  refatorar o que não está quebrado, manter o estilo existente. Código morto não
  relacionado → mencionar, não deletar. Remover imports/variáveis/funções que as próprias
  mudanças tornaram obsoletos.
- Python 3.11; Biopython, pandas, pyyaml. Caminhos **absolutos** no `config.yaml`.
  Threads configuráveis por etapa.
- Se houver abordagem melhor (ferramenta, parâmetro, método estatístico/computacional),
  dizer — discordar quando fizer sentido. Marcar recomendações com "(recomendado)".

## Fluxo de trabalho

- Trabalhar sempre em branch (`feature/...`), nunca direto no `main`.
- Tarefa multi-etapa → declarar plano curto com **checagem verificável por etapa**. Seguir
  sozinho só quando o critério for forte (ex.: "0 stops internos"); critério fraco → parar
  e validar com o usuário.
- Decisões **reversíveis e de baixo risco** (nome de variável, libs equivalentes, formato
  de gráfico): decidir e informar.
- Decisões **irreversíveis ou de alto impacto** (método que muda interpretação biológica,
  exclusão de dados, escolha que define o resto do pipeline): **parar e perguntar antes**.
- Pedir aprovação antes de comandos destrutivos (`git push --force`, `rm`, fechar/mergear
  PR). Mostrar o comando antes de executar qualquer ação de `git`/`gh`.

## Verificação / regressão (caso M. capixaba)

A partir do contig pré-correção, o pipeline deve reproduzir:
ND1 = 882 nt / 294 aa / 0 stops; CYTB = 1143 nt / 381 aa / 0 stops;
ND4 = 1305 nt / 435 aa / 0 stops; ND2 recuperado (≈982 nt).
Mitogenoma final: 19.526 bp, circular, 13 CDS sem stops internos, 2 rRNA, 22 tRNA (MiTFi).
Alvo: `Melipona_capixaba_mitogenome_curated_final_circular_control_region.gb`.

## Ferramentas externas

- **MitoFinder** pode exigir Python 2.7 → tratar como ferramenta externa, com modos
  `python_interpreter` / `conda_env` / `wrapper` no `config.yaml`. Não assumir um env
  `mitofinder_py2`.
- **ARWEN** é opcional (ausente no bioconda do ambiente). Para tRNA, preferir **MiTFi**
  (no caso real, ARWEN superprediz tRNAs).
- **minimap2** por tecnologia: HiFi `-x map-hifi`; ONT `-x map-ont`; PacBio CLR
  `-x map-pb`; Illumina `-ax sr`.

## Editor

vi.

# MitoCurator v0.1-dev

Pipeline diagnóstico para curadoria de mitogenomas montados por long reads, short reads ou dados híbridos.

Esta versão inicial implementa:

- leitura de FASTA ou GenBank;
- checagem de ferramentas externas;
- diagnóstico gene-a-gene de CDS, rRNA, tRNA e regiões intergênicas;
- detecção de stops internos com código genético definido pelo usuário;
- rotação de genomas circulares para iniciar em um gene definido pelo usuário;
- geração de relatórios TSV/Markdown;
- modo `run` all-in-one inicial.

## Instalação do ambiente principal

```bash
mamba env create -f environment.yml
conda activate mitocurator
```

## Uso recomendado

Copie e edite o arquivo:

```bash
cp config.example.yaml config.yaml
```

Use sempre caminhos completos/absolutos no `config.yaml`.

Rodar checagem das ferramentas:

```bash
python -m mitocurator.cli check-tools --config config.yaml
```

Rodar diagnóstico completo inicial:

```bash
python -m mitocurator.cli run --config config.yaml
```

## Curadoria em dois passos: suggest → apply

A etapa de consenso local (`local_consensus`) opera em três modos configuráveis em
`config.yaml`:

| `mode`      | O que faz                                                                 |
|-------------|---------------------------------------------------------------------------|
| `diagnose`  | Só grava o audit log; nenhum arquivo de candidato é escrito.              |
| `suggest`   | Grava candidatos em FASTA + audit log; **sequência do GenBank inalterada**. |
| `apply`     | Aplica candidatos aceitos ao GenBank in-place; escreve `repaired.gb`.     |

**Fluxo recomendado:**

```bash
# 1. Primeira execução — revisar candidatos sem alterar nada
#    No config.yaml: local_consensus.mode: suggest
python -m mitocurator.cli run --config config.yaml

# Revisar:
#   06_local_consensus/audit_log.jsonl  (decisão por gene)
#   06_local_consensus/summary.tsv      (resumo tabular)
#   06_local_consensus/<GENE>/*_candidate.fa  (sequências propostas)

# 2. Segunda execução — aplicar correções aprovadas
#    No config.yaml: local_consensus.mode: apply
python -m mitocurator.cli run --config config.yaml
# Saída: 06_local_consensus/repaired.gb → usado por todas as etapas seguintes
```

Em `mode: suggest` o arquivo `repaired.gb` **não é gerado**; a rotação e o diagnóstico
finais operam sobre o `refined.gb` (sem correções), o que é intencional — o resultado do
diagnóstico reflete o estado pré-correção para facilitar a revisão dos candidatos.

## Estrutura dos diretórios de saída

Os prefixos numéricos dos diretórios de saída (`04_`, `05_`, `06_`, `07_`) identificam
etapas do pipeline, **mas não refletem a ordem de execução**. A sequência real em
`mitocurator run` é:

```
05_refinement  →  06_local_consensus  →  04_rotation  →  07_gene_qc
```

Os números são configuráveis via `config.output.step_dirs` no `config.yaml`. Os valores
padrão foram definidos com lacunas intencionais para acomodar etapas futuras entre eles.

## Observação sobre MitoFinder

O MitoCurator é Python 3. O MitoFinder pode ser chamado separadamente usando:

1. um interpretador Python 2/2.7 explícito;
2. um ambiente conda;
3. um wrapper script.

Exemplo no `config.yaml`:

```yaml
tools:
  mitofinder:
    enabled: true
    mode: python_interpreter
    python2: /usr/bin/python2.7
    script: /home/user/bin/MitoFinder/mitofinder
```

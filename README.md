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

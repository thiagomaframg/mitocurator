# ROADMAP — Auditoria funcional e lógica evidence-driven (branch `codex/01-roadmap-and-functional-audit`)

## Objetivo desta rodada
Documentar o estado funcional atual do MitoCurator e explicitar a lógica de decisão da pipeline com foco em curadoria manual de mitogenoma (caso *Melipona capixaba*), **sem alterar código Python**.

## Escopo confirmado
- Auditoria dos módulos atualmente implementados no repositório.
- Mapeamento do fluxo de decisão operacional já existente.
- Formalização do encadeamento: **problema observado → evidência → diagnóstico → proposta de correção → validação → relatório**.

## Entregáveis desta rodada
1. `docs/current_functional_modules.md`
2. `docs/pipeline_decision_model.md`
3. `docs/mcapixaba_curation_case.md`
4. Este `ROADMAP.md`

## Estado funcional atual (síntese)
- O pipeline funcional está centrado em:
  - `check-tools`
  - `diagnose`
  - `refinement`
  - `rotate`
  - `run` (orquestração inicial)
- Módulos citados no plano amplo (como `read_support.py`, `targeted_consensus.py`, etc.) **não estão presentes neste snapshot** e foram tratados como lacunas de implementação, não como etapas inventadas.

## Prioridades imediatas da documentação
1. Registrar capacidades reais versus capacidades planejadas.
2. Descrever regras de decisão já codificadas (ex.: `internal_stop`, múltiplo de 3, lacunas intergênicas, candidatos de refinamento).
3. Traduzir essas regras para um protocolo de curadoria manual reprodutível do caso *M. capixaba*.

## Critérios de sucesso desta auditoria
- Cobertura dos módulos ativos e suas entradas/saídas.
- Rastreabilidade entre evidência e decisão curatorial.
- Clareza sobre o que é automatizado hoje e o que depende de validação manual.
- Nenhuma alteração em `.py`, `config.example.yaml` e `config.teste.yaml`.

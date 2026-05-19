# Run workflow audit (May 2026)

## 1) Fluxo antigo (mais pesado)

O `run` histórico priorizava um bloco de reconstrução aprofundada, baseado em:

- `read_mapping`
- `variant_evidence`
- `read_support`
- `targeted_extraction`
- `reconstruction_pools`
- `targeted_consensus`
- `candidate_assembly`
- `integrated_report`

Esse encadeamento é útil para análise avançada, porém custa mais CPU/armazenamento e pode ser excessivo para cenários focados em recuperação de gene ausente.

## 2) Módulos chamados atualmente pelo `run` (antes da reorganização)

- check-tools
- anotação inicial (MitoFinder, quando input FASTA)
- refinement
- rotation
- gene QC + annotation assessment
- read mapping + variant evidence
- bloco avançado antigo (condicional por `enabled`):
  - read_support
  - targeted_extraction
  - reconstruction_pools
  - targeted_consensus
  - candidate_assembly
- integrated_report

## 3) Módulos implementados que estavam fora do fluxo principal

- `missing-gene-recovery`
- `missing-gene-assembly-assessment`
- `recovered-contig-annotation`
- `apply-curation` (como subcomando explícito)

## 4) Novo fluxo principal recomendado

Fluxo padrão reorganizado para ser mais enxuto e priorizar recuperação direta de genes ausentes:

0. tool check
1. initial annotation / MitoFinder
2. annotation refinement
3. rotation
4. gene QC + annotation assessment
5. read mapping + read evidence summary (inclui variant evidence quando habilitado)
6. missing-gene-recovery (quando há evidência de CDS ausente/incompleto e entradas de leitura compatíveis)
7. missing-gene-assembly-assessment (quando há assembly local disponível)
8. recovered-contig-annotation (quando há `selected_recovery_contig.fasta`)
9. apply evidence-backed curation (opcional/conservador)
10. final molecule preparation (somente com molécula curada/final disponível)
11. integrated report

### Comportamento de tolerância a falhas esperado no principal

- Se Flye estiver desativado ou indisponível, o `run` registra claramente e segue sem quebrar.
- Se não existir assembly local, o `run` pula assessment/annotation de recovery e continua.
- `apply-curation` permanece opcional por padrão (apenas quando seguro por evidência suficiente).
- `final-molecule-preparation` não é automática sobre contig recuperado.

## 5) Módulos avançados/opcionais

Mantidos como avançados/compatibilidade (não removidos):

- `targeted-extraction`
- `reconstruction-pools`
- `targeted-consensus`
- `candidate-assembly`

Podem ser ativados por configuração quando o estudo exigir reconstrução aprofundada multi-pool/multi-readset.

## 6) Plano de transição

1. Tornar recuperação direta de gene ausente o caminho padrão condicional no `run`.
2. Manter bloco antigo por flags `enabled`, documentado como avançado.
3. Consolidar critérios de aplicação automática de curadoria apenas com suporte robusto.
4. Integrar preparação de molécula final somente quando existir produto curado/final explícito.
5. Atualizar documentação de configuração para formato preferido de `output.outdir`.

## 7) Justificativa: `targeted_consensus` como avançado/opcional

- É um módulo poderoso para cenários complexos, mas computacionalmente mais pesado.
- Em recuperação de gene ausente pontual (como ND2), a trilha direta de recovery + assessment + annotation tende a ser mais eficiente e biologicamente objetiva.
- Mantê-lo opcional evita custo desnecessário no caso padrão sem perder capacidade avançada.

## 8) Justificativa: recovered contig NÃO é automaticamente molécula final

- Um contig recuperado pode resolver genes ausentes e permitir validação funcional local.
- Isso não garante, por si só, que represente a molécula mitocondrial final curada completa (topologia, circularização, consistência global e integração final ainda precisam de validação).
- Portanto, o contig recuperado deve ser tratado como evidência de recuperação/anotação, não como substituto automático da molécula final.

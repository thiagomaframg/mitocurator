# Caso de curadoria — *Melipona capixaba* (modelo operacional)

## Objetivo do caso
Aplicar a lógica evidence-driven do MitoCurator para curadoria manual do mitogenoma de *Melipona capixaba*, mantendo rastreabilidade completa das decisões.

## Fluxo prático no caso

### Etapa A — Problema observado
Entradas a inspecionar:
- `07_gene_qc/problematic_features.tsv`
- `05_refinement/expected_gene_set.tsv`
- `05_refinement/missing_gene_candidates.tsv`
- `05_refinement/cds_refinement_candidates.tsv`
- `07_gene_qc/intergenic_regions.tsv`

Sinais de alerta esperados no caso:
- CDS com stop interno.
- CDS fora de múltiplo de 3.
- gene mitocondrial esperado ausente/duplicado.
- lacuna intergênica longa/AT-rich compatível com região de controle ou com gene não anotado.

### Etapa B — Evidência
Para cada problema, registrar no caderno de curadoria:
- gene/região;
- coordenadas atuais;
- tipo de evidência (stop interno, frame, ausência, candidato ORF, lacuna);
- prioridade pelo `decision_hint`.

### Etapa C — Diagnóstico
Classificar a hipótese principal:
1. fronteira CDS imprecisa;
2. indel/frameshift potencial;
3. gene ausente na anotação;
4. região não codificante mal anotada;
5. duplicação anotacional potencial.

### Etapa D — Proposta de correção
Ajustes permitidos no processo curatorial (fora desta auditoria documental):
- editar fronteira de CDS com melhor candidato de refinamento;
- incorporar candidato de gene ausente quando houver sustentação;
- manter status “incerto” quando evidência for fraca.

### Etapa E — Validação
Após cada ajuste:
1. rodar novamente `diagnose`;
2. comparar antes/depois em `problematic_features.tsv`;
3. verificar consistência no `diagnostic_summary.md`.

Critério de aceitação local:
- diminuição objetiva de inconsistências sem introduzir novos conflitos óbvios.

### Etapa F — Relatório final
Relatar cada decisão no formato:
- **Problema**
- **Evidência**
- **Diagnóstico**
- **Correção proposta/aplicada**
- **Validação (antes/depois)**
- **Status** (aceita, pendente, rejeitada)

## Limites atuais do caso no snapshot
- Não há, ainda, validação por suporte direto de leituras no MitoCurator atual.
- Não há módulos ativos de extração/consenso direcionado para fechar regiões ambíguas automaticamente.
- Portanto, a curadoria de *M. capixaba* neste estágio depende de revisão manual assistida pelos TSVs diagnósticos e de refinamento.

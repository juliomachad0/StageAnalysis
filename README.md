# Stage Analysis
StageAnalysis

## CLASSE DE AVALIAÇÃO DE DV

De modo a avaliar o dV para diferentes motores e em diferentes configurações, foi desenvolvida
]uma classe chamada deltaV. Essa classe permite a configuração dos estágios de um foguete e o 
cálculo automático do dv associado a tal configuração, contando com a opção de considerar massas
que serão alijadas em determinados estágios, como, por exemplo, a coifa que é descartada 
geralmente na queima do penúltimo/antepenúltimo estagio.

### métodos:
#### add_stage:
adiciona estagio - configuração (isp, mp, me, efic) e posição (de baixo para cima).
Onde: mp-massa de propelente, me-massa da estrutura do estágio (dry mass), efic-eficiência máxima.

#### remove_stage:
remove estagio, informa-se a posição do mesmo.

#### update_stage: 
atualiza estagio. Informa-se a nova configuração de (isp, mp, me, efic) e 
posição do estágio. O método utiliza os métodos remove_stage e add_stage

#### update_jettison_mass: 
atualiza a massa de alijamento (massa, posição do estágio onde é liberada)

#### update_payload_mass:
atualiza a carga paga

#### show_up_stages: 
mostra a configuração de estágios posicionados ordenadamente (mp, me, isp)

#### deltaV: 
calcula o deltaV do estágio da configuração de estágios informada.

#### get_optimum_prop_qtd: 
calcula a quantidade de propelente necessário a um dado estágio para que
o foguete atinja um DeltV fornecido como alvo.

#### get_max_payload_mass: 
Quando o DeltaV da configuração de estágios é maior que o DeltaV 
fornecido como referência/alvo, o método aumenta a carga paga até que DeltaV da
configuração se iguale ou fique menor que o DeltaV alvo.

### Principais parâmetros:

#### stages: 
lista dos estágios, ex: [[m prop, m estrutura, isp],[mesma coisa pro segundo estagio] ...]

#### jettison_mass: 
massa que vão ser alijadas (ex coifa) na forma 
[numero do estágio a ser alijada, massa a ser alijada], ex: [[2, massa1],[3,massa2], 
ex: massa coifa a ser alijada no segundo estágio: jettison_mass = [2, mcoifa]. 
Se for para ser considerada até o final (ex: baterias) colocar 'n' 
no primeiro parâmetro da lista, ex(massa baterias): ['n', mbaterias]

payload_mass: autoexplicativo, isp: autoexplicativo, deltav_reference/deltaV_target: autoexplicativo, 
position/stage_position: posição do estágio considerado, 
de baixo pra cima (exemplo número do estágios do vlm: 
[S50, S50, terceiro liquido] -> númeração: [1,2,3], 1: S50, 2: S50, 3: terceiro liquido

mp: massa do propelente, me: massa das estruturas (calculada a partir das 
eficiências mássicas ou servindo de base para o cálculo de tais eficiências. 
Massa dos tanques ou qualquer outra que varie com a quantidade de propelente),
efic: eficiência mássica (mp/(mp+me)),

## CLASSE DE EFICIENCIA DE STAGIO



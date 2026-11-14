# Especificação da física do notebook

Estado: versão 8 da especificação (seção 15, E8, α por faixas estreitas de Teff, aprovada e implementada). Versão 7 da especificação (E7c: CR-01, análise principal com energia escalada pela luminosidade, seção 14). Versão 6 da especificação (E7b: tamanho mínimo da cauda na varredura de xmin e achados da prévia). Versão 5 da especificação (correções após a E7a: nomes das tabelas do VizieR, origem de Teff e log g do Yang, procedência BTSL). Versão 4 da especificação (E6 fechada na versão 3, pendências conferidas nos PDFs na abertura da E7). Versão 3 da etapa E6 (rascunho 2 aprovado pelo usuário, com idioma inglês para as saídas e o notebook). As propostas da seção 11 foram aprovadas pelo usuário, com RAD_MAX definido em relação à mediana (seção 3.2). Estrelas F entram como extensão secundária (seção 13), decisão do usuário. Este documento define o que o notebook calcula e por quê. Não contém código de análise. Fatos da literatura vêm com chave bibliográfica e página do PDF. Itens marcados como "escolha de projeto" não têm fonte na literatura lida e são decisões de desenho a aprovar.

Convenções do texto: sem travessões, sem ponto e vírgula, intervalos escritos com "to", energias em erg, expoente α sempre positivo, com dN/dE proporcional a E^(-α).

## 1. Pergunta e hipóteses

Pergunta: o expoente α da distribuição de energia de flares estelares é o mesmo para anãs M, K e G?

Hipóteses principais, para as classes c em {M, K, G}:

* H0: α_M = α_K = α_G (expoente comum).
* H1: pelo menos um dos três expoentes difere dos demais.

Hipóteses secundárias, uma por classe:

* H0c: α_c = 3/2, valor previsto pelo modelo fractal-difusivo de criticalidade auto-organizada (FD-SOC) para a energia, com d = 3, β = 1 e D_3 = 2 (aschwanden2014 p8, Eqs. 25 a 27, aschwanden2016 p18, Eqs. 22 a 24).
* H1c: α_c ≠ 3/2.

Hipóteses da extensão com estrelas F (secundárias, declaradas de antemão, seção 13):

* H0F: α_F = α_G = α_K = α_M. H1F: pelo menos um difere.
* H0FG: α_F = α_G (comparação por par).

Limites do enquadramento (escolhas congeladas):

* O notebook não afirma pertencimento a uma classe de universalidade de sandpile. Os expoentes do modelo BTW em 2D não são conhecidos analiticamente e a interpretação é controversa (lubeck2004 p113). Os valores de modelos de sandpile em 2D compilados vão de 0,98 a 1,7 (aschwanden2016 p28 e p29, Tabela 1).
* Os valores "1,8 para BTW 2D" e "2 para BTW" do documento original não são usados. O 1,8 vem da inclinação do fluxo de pico de flares solares em raios X duros (lu1991 p1 e p2, crosby1993 p282).
* Uma lei de potência aparente não prova criticalidade em sistemas sem conservação (bonachela2009 p1), e a evidência observacional de SOC é difícil de obter (watkins2016 p27). A discussão do manuscrito deve respeitar isso.

## 2. Dados

### 2.1 Fontes e acesso

O notebook baixa os catálogos do VizieR por código (astroquery) e não versiona dados brutos. A cada execução registra versão do catálogo, data de acesso, número de linhas e um checksum do conteúdo baixado.

| Papel | Catálogo VizieR | Tabela | Conteúdo usado |
|---|---|---|---|
| Flares | J/ApJS/241/29 (Yang e Liu 2019) | table2 | KIC, trimestre Q, Begin, End, logE (162262 linhas) |
| Estrelas com flares | J/ApJS/241/29 | table1 | KIC, Nf, Teff e logg do KIC, Prot, NS (3420 linhas) |
| Parâmetros estelares | J/ApJS/229/30 (Mathur et al. 2017) | catalog (no VizieR as table3 e table4 do ReadMe aparecem fundidas nesta tabela) | KIC, Teff, log(g), Rad, Mod e intervalos de confiança de 1 sigma E_ e e_ (197096 linhas) |

Notas sobre as fontes:

* A table1 do Yang reflete o erratum de 2020. O erratum corrige apenas a razão L_flare/L_bol e as Eqs. 1 e 2 do artigo (yang2020erratum p1 e p2). Teff, logg, Nf, Prot e logE não mudam. O notebook não usa L_flare/L_bol.
* O Yang usa fluxo SAP de cadência longa, trimestres Q1 a Q17 (yang2019 p3).
* O Yang cobre estrelas de 2500 K a 11287 K, com logg de 2,2 a 5,4 (ReadMe), isto é, inclui subgigantes e gigantes.
* A nota da Tabela 1 do Yang diz que Teff e log g são dados pelo KIC (yang2019 p5), e a p9 diz que os parâmetros são os do DR25, calibrados por isócronas de Dartmouth (Mathur 2017). Os dados confirmam a segunda leitura. Nas 3302 estrelas casadas, Teff e log g da table1 coincidem com os do Mathur (mediana da diferença 0 K e 0,000 dex, desvio padrão de 13,5 K e 0,004 dex, 3 estrelas com diferença maior que 50 K). O catálogo antigo do KIC subestima Teff em até 200 K e superestima log g de estrelas tipo solar (mathur2017 p2), de modo que usar os valores DR25 evita esse viés. A table1 não traz raio, intervalos de confiança nem procedência. O Mathur fornece tudo isso, e por isso o cruzamento continua necessário, embora Teff e log g já coincidam.

### 2.2 Verificações de integridade (executadas e relatadas)

1. Número de linhas: 3420 na table1, 162262 na table2, 197096 na table4 do Mathur.
2. Faixa de logE entre 29,3 e 38,6 (ReadMe) e ausência de valores não finitos.
3. Toda KIC da table2 existe na table1, e Nf da table1 igual ao número de flares por KIC na table2.
4. Fração das 3420 estrelas com contraparte no Mathur. Relatar as que faltam e excluí-las com contagem explícita.
5. Comparação Teff do Yang (KIC) contra Teff do Mathur: mediana e dispersão da diferença, por classe.
6. Flares duplicados (mesma KIC e mesmo Begin): contar e relatar, sem remover sem decisão.

### 2.3 Energia

A energia de cada flare é E = 10^logE erg, exatamente como publicada. O cálculo do Yang assume corpo negro de 9000 K, com área do flare obtida da amplitude, raio estelar, Teff e resposta do Kepler (yang2019 p4 e p5). O erro de energia é de cerca de 60 por cento e a cadência longa subestima energias em cerca de 25 por cento (yang2019 p5, yang2018 p1).

Propriedade usada na análise: uma constante multiplicativa em todas as energias não altera α, porque desloca a distribuição na escala logarítmica sem mudar sua inclinação. Um fator que dependa da classe poderia deslocar a faixa de energia comparada, e por isso entra como sistemática (seção 7).

### 2.4 Efeitos de seleção do catálogo que o notebook deve declarar

* Sem teste de injeção e recuperação e sem limite de completude por estrela (yang2019 p5). A completude é tratada por xmin (seção 5.3).
* Estrela com vários trimestres só tem flares aceitos se aparecerem em pelo menos 2 trimestres, e estrela com um único candidato só entra se E for maior que 1e34 erg e duração de pelo menos 2 h (yang2019 p3 e p4). Esse critério depende da energia e atinge o regime de baixa energia dessas estrelas. A sensibilidade da seção 7 exclui estrelas com Nf igual a 1.
* Validação visual com escolhas arbitrárias para candidatos marginais (yang2019 p4).
* Estrelas com vizinha a menos de 12 segundos de arco (152, código NS igual a p) estão sinalizadas e não removidas (yang2019 p4, ReadMe). Há sensibilidade que as exclui.
* A cadência longa perde cerca de 60 por cento dos flares detectados na cadência curta (yang2018 p1), o que afeta a região de baixa energia.

## 3. Amostra e filtro de anãs

### 3.1 Funil de seleção

O notebook imprime uma tabela com o número de estrelas e de flares após cada passo:

1. Todos os flares da table2 (162262).
2. Flares de estrelas com contraparte no Mathur.
3. Flares de estrelas que passam o filtro de anãs (seção 3.2).
4. Flares de estrelas com Teff dentro de uma das três classes (seção 4).
5. Flares com logE finito, que formam a amostra analisada. A implementação inclui o passo 3b, que remove estrelas pelas opções de exclusão das sensibilidades e não faz nada na configuração principal.

### 3.2 Filtro de anãs

Motivação: cerca de 40 por cento das estrelas de superflares dos primeiros estudos eram subgigantes ou gigantes (notsu2019 p27, okamoto2021 p9 e p10), e o DR25 inclui 43,5 por cento mais subgigantes que o catálogo anterior (mathur2017 p1). Dos alvos Kepler com raios da Gaia DR2, cerca de 67 por cento são de sequência principal, 21 por cento subgigantes e 12 por cento gigantes vermelhas (berger2018 p1).

Papel do raio nesta especificação: o raio serve só para selecionar anãs, isto é, para excluir subgigantes e gigantes da amostra de estrelas com flares. As estrelas com flares e as energias vêm do Yang, e o Mathur entra apenas pelo cruzamento por KIC. O raio do Mathur não é usado para recalcular energias nem aparece em nenhuma outra análise. Observação: a energia do Yang depende do raio estelar (yang2019 p4) e as luminosidades do artigo usam raios do KIC (yang2019 p9), então uma estrela cujo raio do KIC difere do raio do Mathur carrega erro extra na energia. O filtro de anãs reduz esse risco porque remove justamente as estrelas evoluídas mais mal classificadas pelo KIC. A fonte do raio usada no cálculo da energia não está explícita nas páginas lidas.

Critério (decisão DEC-17): usar log g e raio do Mathur (colunas log(g) e Rad da table4). Parâmetros:

* LOGG_MIN. Valor proposto: 4,0. Fonte do valor: definição de anãs G com log g maior que 4,0 (shibayama2013 p6 e p7, notsu2013 p3 e p5). Variante: 4,2 (walkowicz2011 p4).
* RAD_MAX. Nenhum PDF da coleção fornece um limite numérico de raio para anãs. Os artigos de Mathur, Huber e Berger descrevem classes evolutivas sem fixar um corte simples (berger2018 p7 usa trilhas evolutivas no plano Teff e raio). Definição adotada (escolha de projeto aprovada): o corte é relativo à mediana. Para cada intervalo de Teff de 200 K (escala da incerteza do KIC, brown2011 p2), calcula-se a mediana R_med do raio (coluna Rad do Mathur) das estrelas com flares ligadas ao Mathur, antes do filtro. A estrela é mantida se Rad <= RAD_FATOR × R_med do seu intervalo. Intervalos com menos de 30 estrelas são fundidos ao vizinho. A mediana representa a sequência principal porque cerca de 67 por cento dos alvos Kepler são de sequência principal (berger2018 p1). RAD_FATOR proposto: 1,5, com variantes 1,25 e 2,0 (o fator não tem fonte e é escolha de projeto).
* Estrelas sem classificação: o Yang diz que estrelas não classificadas têm parâmetros iguais aos do Sol e ficam fora dos seus estudos (yang2019 p9, citando Mathur). No Mathur, porém, a table4 traz só as 197096 estrelas com parâmetros derivados, e as 2942 estrelas não classificadas ficam fora dela, sem Teff, log g nem metalicidade (mathur2017 p7). A atribuição de valores solares a não classificadas pertence ao catálogo anterior (huber2014 p13, procedência Solar, 7 por cento da amostra daquele catálogo). Regra: flare sem contraparte na table4 é excluído e contado no passo 2 do funil. No conjunto real, 118 das 3420 estrelas do Yang não têm contraparte e saem nesse passo. Procurar triplas repetidas de Teff, log g e Rad não serve como verificação, porque o ajuste de isócronas gera muitos valores idênticos que não são valores solares padrão. A coluna Mod do Mathur indica a procedência: DSEP para ajuste de isócronas, BTSL para 235 estrelas de anãs frias ou fora da grade, com valores obtidos por ajustes polinomiais a modelos BT-Settl e sem distância nem extinção (mathur2017 p6), e MULT para um punhado de estrelas (3 no catálogo). Entre as estrelas com flares, 54 são BTSL. O notebook as conta e as mantém, e uma sensibilidade as exclui (S7).
* Incertezas: o Mathur fornece intervalos de 1 sigma para log g e raio. O filtro usa o valor central. Uma variante exige que o corte seja satisfeito também no extremo do intervalo.

O notebook faz a análise completa em três configurações: sem filtro, filtro apenas por log g e filtro por log g e raio. O resultado principal usa a configuração aprovada na seção 11. As outras entram como sensibilidade.

## 4. Classes de Teff

Teff vem do Mathur. Fronteiras primárias pela escala de anãs de Pecaut e Mamajek (pecaut2013 p9, Tabela 5): G0V em 5920 K, K0V em 5280 K, M0V em 3850 K, M4V em 3200 K. Arredondamento adotado (decisão D-A): F e G em 5900 K, G e K em 5300 K, K e M em 3850 K.

| Classe | Intervalo de Teff |
|---|---|
| G | 5300 K <= Teff < 5900 K |
| K | 3850 K <= Teff < 5300 K |
| M | Teff < 3850 K |

Estrelas mais quentes que 5900 K não entram na análise principal. A fronteira com F delimita G e define a classe F da extensão (seção 13).

Variantes de fronteira (seção 7):

* Fronteiras do Yang: 6000 K, 5000 K e 4000 K (yang2019 p7, Tabela 3).
* Deslocamento simultâneo de todas as fronteiras por menos 200 K e mais 200 K, pela incerteza do KIC e das escalas de Teff (brown2011 p2, huber2014 p2).
* Valores de Pecaut sem arredondamento: 5920 K, 5280 K e 3850 K.
* Subdivisão da classe M em M3200 (3200 K <= Teff < 3850 K) e M totalmente convectiva (Teff < 3200 K), seguindo a transição em torno de M4 (hawley2014 p2, yang2019 Fig. 3, gunther2020 p12). Os 3200 K coincidem com M4V em pecaut2013 p9.

## 5. Estimadores e testes por classe

Para cada classe, com energias x_i e n_tail = número de flares com x >= xmin.

### 5.1 Modelo e estimador

Lei de potência contínua acima de xmin:

p(x) = ((α menos 1) / xmin) (x / xmin)^(-α), com α > 1.

Estimador de máxima verossimilhança (clauset2009 p5, Eq. 3.1):

α_hat = 1 + n_tail / soma de ln(x_i / xmin).

Erro padrão de referência (clauset2009 p5, Eq. 3.2): σ = (α_hat menos 1) / raiz(n_tail). Ele serve apenas como comparação. As incertezas reportadas vêm do bootstrap (seção 5.5).

A implementação usa o pacote powerlaw (alstott2014 p3), em modo contínuo (alstott2014 p7). Os valores de logE têm 4 casas decimais e a energia é tratada como contínua. O notebook verifica quantos valores repetidos existem e relata.

Regressão por mínimos quadrados em histograma log-log não é usada como estimador principal, porque produz erros sistemáticos (clauset2009 p5). Aparece só na reprodução do Yang (seção 7, S8).

### 5.2 Escolha de xmin

xmin é o valor que minimiza a distância KS entre os dados com x >= xmin e o ajuste (clauset2009 p11, Eq. 3.9, alstott2014 p5). A varredura percorre os valores únicos dos dados. Com cerca de 1e5 flares por classe a varredura completa é cara. Escolha de projeto: o notebook usa uma grade de candidatos a xmin em passos de 0,02 dex, e verifica em subamostras que a grade e a varredura completa dão o mesmo α dentro da incerteza de Hill (verificado na E7b, diferença máxima de 0,008 em dados sintéticos e 0,007 nos dados reais). A varredura exige pelo menos 1000 pontos na cauda, o valor a partir do qual Clauset relata bons resultados (clauset2009 p13). Sem essa exigência, a distância KS mínima do K e do G caía em xmin muito altos, com 608 e 834 pontos na cauda, e o critério de validade da seção 8.1 falhava por construção. Nas validações com amostras sintéticas pequenas o mínimo é de 100 pontos.

KS é preferido ao BIC, que tende a subestimar xmin (clauset2009 p11), e a estatística de Anderson-Darling dá xmin grande demais (clauset2009 p13). Aschwanden 2015 argumenta que KS tende a dar limiar alto e usa minimização de qui-quadrado em histogramas (aschwanden2015 p8). Esse método não é MLE e não entra na análise principal. Pode aparecer como comparação.

### 5.3 Completude e xmin

Sem injeção e recuperação no catálogo, xmin do KS serve de proxy de completude. Sensibilidade: xmin fixo em 1e33, 1e33,5 e 1e34 erg para todas as classes. A ausência de completude medida é declarada como limitação.

### 5.4 Truncamento

Efeitos que afastam a distribuição da lei de potência ideal: limiar físico, amostragem incompleta de eventos pequenos, contaminação de fundo e truncamento por tamanho finito (aschwanden2015 p1). O rollover por subamostragem achata a parte baixa e o tamanho finito aumenta a inclinação no alto, com quase nenhuma faixa intermediária confiável (aschwanden2015 p13). Em amostras com tempo finito, a distribuição cumulativa de uma lei de potência finita cai como exponencial perto do maior evento (aschwanden2016 p82 e p83).

Modelos ajustados acima de xmin:

* M1, lei de potência pura (xmax infinito).
* M2, lei de potência com truncamento superior duro em xmax. A densidade normalizada é p(x) = (α menos 1) x^(-α) / (xmin^(1-α) menos xmax^(1-α)), para xmin <= x <= xmax. O pacote aceita xmax (alstott2014 p5 e p6). Valor de xmax: a maior energia observada na classe.
* M3, lei de potência com corte exponencial, para comparação.
* M4, lognormal, para comparação.

Resultado principal: M2 (decisão D-D inclui truncamento superior). M1 é reportado lado a lado. Escolha de projeto: confirmar M2 como principal na seção 11.

Sensibilidade ao truncamento: xmax igual ao quantil 99,9 por cento e 99 por cento das energias da classe, para checar a influência dos maiores eventos. As energias máximas podem ter erro (davenport2016 p13 alerta para energias de 1e39 erg por erro de de-trending ou de luminosidade).

### 5.5 Incertezas por bootstrap

Bootstrap não paramétrico, com reajuste completo (xmin incluído) em cada réplica. Número de réplicas: B_BOOT (padrão 1000).

* Bootstrap por estrela (principal): reamostra estrelas com reposição dentro da classe e leva todos os flares de cada estrela sorteada. Escolha de projeto: flares da mesma estrela não são independentes, então o reamostrador por estrela é o mais honesto. Quando a estrela muda de classe por incerteza em Teff, isso é tratado na sensibilidade de fronteiras, não aqui.
* Bootstrap por flare (comparação): reamostra flares com reposição.
* Perturbação de energia: em cada réplica, cada logE recebe ruído gaussiano de desvio padrão 0,2 dex (cerca de 60 por cento em energia, yang2019 p5, decisão D-C). Executar com e sem a perturbação e relatar a diferença.

Saídas: mediana, intervalo de percentis 16 a 84 e intervalo de 99 por cento de α e de xmin, por classe.

### 5.6 Aderência

Teste de aderência semiparamétrico por Monte Carlo (clauset2009 p15 a p17): ajusta, calcula KS, gera conjuntos sintéticos com mesmos α e xmin (abaixo de xmin usa os dados empíricos), reajusta cada sintético incluindo xmin e conta a fração de KS sintéticos maiores que o empírico. Critério: p > 0,1 sustenta a lei de potência, p pequeno a descarta (clauset2009 p17 e p18). Número de sintéticos: N_GOF (padrão 200). Para que o p-valor tenha precisão de cerca de ε, são necessários ao menos 1/4 ε^(-2) conjuntos (clauset2009 p17), então N_GOF igual a 200 dá ε de cerca de 0,035, suficiente para a decisão em torno de 0,1 e não para citar p com duas casas. Clauset usou de 1000 a 10000 conjuntos (clauset2009 p17). O notebook relata ε e, no modo completo, pode elevar N_GOF. Para n_tail pequeno (abaixo de cerca de 100), p alto exige cautela (clauset2009 p18).

Observação de projeto: com n_tail da ordem de 1e4 ou mais, o teste pode rejeitar desvios fisicamente pequenos. Por isso o critério de decisão da seção 8 nunca depende apenas do p de aderência. O notebook relata também a diferença máxima KS D.

### 5.7 Comparação com alternativas

Razão de verossimilhanças R de M1 (ou M2) contra M3 e M4, com o mesmo xmin nos dois modelos. O sinal de R é confiável quando o p do teste de Vuong é menor que 0,1 (clauset2009 p19 e p20). Para modelos aninhados (lei de potência contra lei de potência com corte exponencial) é preciso o teste modificado (clauset2009 p19 e p20, Apêndice C). O pacote powerlaw calcula R e p com distribution_compare (alstott2014 p3 e p8), aceita o argumento nested e detecta sozinho quando um modelo é a versão aninhada do outro, aplicando a correção do p-valor (alstott2014 p12). Na comparação lei de potência contra lei de potência truncada, o pacote chama de truncamento o corte exponencial (alstott2014 p9), que é o M3. O argumento normalized_ratio normaliza R pelo desvio padrão (alstott2014 p8).

## 6. Teste de igualdade de α entre M, K e G

### 6.1 Estatística

Para cada classe c há n_c flares com x >= xmin_c. Com xmin_c fixos nos valores estimados, a verossimilhança de cada classe é a do modelo M2 (ou M1). Defina:

* l_sep, soma das log-verossimilhanças máximas com α_c livres.
* l_com, log-verossimilhança máxima com um único α comum às três classes.
* Λ = 2 (l_sep menos l_com), estatística de razão de verossimilhanças para modelos aninhados.

Para M1, o α comum tem forma fechada: α_com = 1 + N / soma_c S_c, com N a soma dos n_c e S_c a soma de ln(x_i / xmin_c) na classe c. Para M2 a maximização é numérica.

Os testes por pares M contra K, K contra G e M contra G usam a mesma construção com duas classes.

### 6.2 Calibração por simulação

Escolha de projeto, porque xmin é estimado dos dados e a distribuição de Λ sob H0 não é um qui-quadrado simples:

1. Ajustar o modelo com α comum às três classes.
2. Para cada réplica sob H0 (B_NULL, padrão 500), gerar para cada classe uma amostra semiparamétrica como em clauset2009 p15 a p17: acima de xmin_c, lei de potência com α comum, abaixo, reamostra dos dados empíricos, com o mesmo n_c.
3. Reajustar cada classe (incluindo xmin) e recalcular Λ.
4. p-valor igual à fração de réplicas com Λ simulado maior ou igual ao observado.

### 6.3 Medidas de tamanho de efeito

* Diferenças Δα_MK, Δα_KG e Δα_MG com intervalos de bootstrap por estrela.
* Valor de α comum e sua incerteza.
* Quando as classes têm xmin muito diferentes, a comparação mistura intervalos de energia diferentes (o Yang usa faixas fixas distintas, 1e33,5 a 1e35 erg para K e M e 1e34,5 a 1e36 erg para G, yang2019 p7 e p8). Por isso há uma análise em janela comum, com xmin igual ao maior xmin_c entre as classes, e as análises com xmin fixo da seção 5.3.

### 6.4 Comparação com a literatura

O notebook coloca os resultados ao lado dos valores publicados, rotulando o método de cada um:

| Fonte | Classe | α | Método |
|---|---|---|---|
| yang2019 p7 e p8, Fig. 3 (valores lidos na figura) | G | 1,96 ± 0,04 | OLS bisector, histograma diferencial |
| | K | 1,78 ± 0,02 | idem |
| | M | 2,13 ± 0,05 | idem |
| | M totalmente convectiva | 2,09 ± 0,10 | idem |
| shibayama2013 p15 | G (todas / rotação lenta) | 2,2 / 2,0 | ver artigo |
| okamoto2021 p16 a p19, Tabela 4 | G 5100 a 5600 K | 1,5 a 2,2 conforme rotação | FFD corrigida por sensibilidade |
| | G 5600 a 6000 K | 1,8 a 2,7 conforme rotação | idem |
| wu2015 p1 | G | 2,04 ± 0,17 | lei de potência, 5e34 a 1e36 erg |
| tu2020 p1 e p2 | tipo solar (TESS) | 2,16 ± 0,10 | lei de potência |
| aschwanden2016 p82 | Kepler, G5 | 2,04 ± 0,13 | ver artigo |
| FD-SOC | | 3/2 | aschwanden2014 p8 |

Os valores do Yang foram conferidos em zoom na Figura 3 (yang2019 p7): A 1,12 ± 0,08, F 2,11 ± 0,09, G 1,96 ± 0,04, K 1,78 ± 0,02, M 2,13 ± 0,05, M totalmente convectiva 2,09 ± 0,10, gigantes 1,90 ± 0,10. A mesma figura traz os números de flares do Yang por tipo: A 583, F 8869, G 55259, K 47112, M 50439, gigantes 6496. Eles servem de expectativa para n_tail, que pode ficar abaixo de 1000 para F depois de xmin e dos filtros. Okamoto reporta o índice como valor negativo na convenção dE/dN proporcional a E^α, então aqui entram em módulo.

## 7. Sensibilidades

Cada variante repete o ajuste, o teste global e os testes por pares. O notebook registra para cada uma α por classe, xmin, n_tail, p global, veredito e contagem de classes com n_tail >= 1000 (clauset2009 p13).

| Código | Variação |
|---|---|
| S0 | Configuração principal |
| S1 | xmin fixo em 1e33, 1e33,5 e 1e34 erg (decisão D-E) |
| S2 | Fronteiras de Teff: Yang, menos 200 K, mais 200 K, Pecaut sem arredondar, M subdividida |
| S3 | Filtro de anãs: sem filtro, só log g (4,0 e 4,2), log g e raio, filtro nos extremos de 1 sigma |
| S4 | Energia: perturbação lognormal de 0,2 dex nas réplicas. Fator multiplicativo de 1,25 aplicado a uma classe de cada vez (e seu inverso), como teste de que o pipeline é invariante a escala quando o xmin é livre |
| S5 | Truncamento: M1, M2, xmax nos quantis 99,9 e 99 por cento |
| S6 | Heterogeneidade por estrela (detalhada abaixo) |
| S7 | Exclusões: estrelas com Nf igual a 1, estrelas com NS igual a p, estrelas sem período de rotação (Prot igual a menos 1), estrelas com Mod diferente de DSEP |
| S8 | Reprodução do Yang: bissetriz de OLS (isobe1990) no histograma diferencial com faixas fixas de energia (1e34,5 a 1e36 erg para G e 1e33,5 a 1e35 erg para K e M, yang2019 p7 e p8) |
| S9 | Janela comum de energia (xmin igual ao maior xmin entre classes) |

Detalhe de S6 (decisão D-F): o resultado principal usa flares agrupados por classe. Como verificação de heterogeneidade, o notebook seleciona estrelas com pelo menos N_MIN flares acima do xmin da classe (N_MIN proposto: 50, escolha de projeto). Estima α por estrela com xmin fixo no da classe e compara a dispersão observada com a esperada pelo erro de amostragem (alpha_hat menos 1) / raiz(n) (clauset2009 p5), com um teste de heterogeneidade do tipo Q de Cochran por classe (escolha de projeto). Davenport usou estrelas com pelo menos 100 candidatos (davenport2016 p7), mas aquele conjunto usa mínimos quadrados e contém não anãs.

Fórmula usada em S8 (isobe1990 p106, Tabela 1): com x igual a logE do centro do bin, y igual a log10 de dN/dE por ano, e as somas Sxx, Sxy e Syy dos desvios em relação às médias, beta1 = Sxy / Sxx, beta2 = Syy / Sxy e a inclinação da bissetriz beta3 = (beta1 + beta2)^(-1) [beta1 beta2 menos 1 + raiz((1 + beta1^2)(1 + beta2^2))]. O Yang usa histograma diferencial com barras de erro de Poisson (yang2019 p7). O valor de α é a inclinação, em módulo.

Detalhe de S8: o objetivo é mostrar quanto do resultado do Yang vem do estimador e das faixas fixas. Não se exige igualdade com os valores publicados. Os artigos do Yang e de outros autores reconhecem que α depende da faixa de energia, do número de estrelas e do instante (yang2019 p8).

## 8. Critérios de decisão

### 8.1 Validade do ajuste por classe

Uma classe é considerada bem caracterizada quando:

1. n_tail >= 1000 (clauset2009 p13).
2. O p de aderência é maior que 0,1 (clauset2009 p17 e p18) ou, se menor, o desvio KS D e a comparação por R (Vuong com p < 0,1) mostram que a lei de potência segue sendo o melhor modelo entre os testados.
3. A faixa dinâmica de energia da cauda (log10 de xmax/xmin) é relatada. Nenhum limite mínimo de faixa tem fonte na literatura lida, então fica como relato e não como critério.

Classes que não cumprirem 1 ficam marcadas e os resultados são lidos com essa ressalva.

### 8.2 Veredito sobre H0

Para a análise principal definida na seção 14 (CR-01) valem as regras daquela seção. O texto abaixo vale para a linha de base em energia como publicada.

Nível de significância: 0,05 para o teste global. Comparações por pares com correção de Holm (escolha de projeto).

* H0 não rejeitada: p global >= 0,05 na configuração principal.
* H0 rejeitada: p global < 0,05 na configuração principal.

Margem de relevância: DELTA_EQ proposto em 0,10 em α (escolha de projeto, a aprovar). Se H0 for rejeitada mas todos os |Δα| forem menores que DELTA_EQ, o relato diz que a diferença é estatisticamente detectável e fisicamente pequena.

### 8.3 Robustez

O veredito da configuração principal é chamado robusto se for o mesmo em todas as variantes de S1 a S5, S7 e S9 que mantenham n_tail >= 1000 nas três classes. Caso contrário, o veredito é descrito como dependente de escolhas de análise, com a lista das variantes que o invertem. S8 e S6 são diagnósticos e não entram na regra.

### 8.4 Comparação com 3/2

Para cada classe, H0c é rejeitada quando 3/2 fica fora do intervalo de 99 por cento do bootstrap por estrela de α. O notebook relata α_c menos 3/2 e o intervalo. Não há afirmação sobre classe de universalidade. As explicações possíveis para α acima de 3/2, já presentes na literatura, são a relação não linear entre energia bolométrica e energia dissipada (aschwanden2016 p83) e os efeitos de truncamento (aschwanden2015 p13). O notebook não decide entre elas.

## 9. Saídas

### 9.1 analysis_results.md

Gerado pelo notebook a cada execução, em inglês, em ordem fixa:

1. Cabeçalho: data, versões de Python, numpy, scipy, powerlaw, astropy, astroquery e matplotlib, sementes, modo de execução (completo ou rápido), checksums dos dados.
2. Integridade dos dados (seção 2.2).
3. Funil de seleção (seção 3.1).
4. Tabela principal por classe: n_estrelas, n_flares, xmin, n_tail, α com incerteza, p de aderência, KS D, R e p de Vuong contra M3 e M4, modelo preferido.
5. Teste global e testes por pares: Λ, p, Δα com intervalos.
6. Comparação com 3/2 por classe.
7. Tabela de sensibilidades S0 a S9 e veredito de robustez.
8. Heterogeneidade por estrela (S6).
9. Tabela de comparação com a literatura (seção 6.4).
10. Limitações automáticas: classes com n_tail abaixo de 1000, variantes que invertem o veredito.

### 9.2 results_macros.tex

Arquivo LaTeX com um \newcommand por número usado no manuscrito. Regras:

* Nomes só com letras (sem dígitos nem sublinhado), no padrão grandeza seguida de sufixo de classe: AlphaM, AlphaK, AlphaG, AlphaErrM, XminM, NtailM, NflaresM, NstarsM, PgofM, DeltaAlphaMK, PglobalTest, NvariantsTotal, NvariantsSameVerdict.
* Formatação: ponto decimal, sem separador de milhar, α e incertezas com 2 casas, p-valores com 3 casas, log10 de xmin com 2 casas, contagens como inteiros.
* Cabeçalho do arquivo com data, modo de execução e checksums. O arquivo é regenerado inteiro, nunca editado à mão.
* O manuscrito usa apenas esses macros para números, de modo que texto e resultado não divirjam.

### 9.3 Figuras

Estilo comum (decisão DEC-18, normas da EPL):

* Layout Super Mongo reproduzido no matplotlib: preto e branco, marcas dos eixos para dentro nos quatro lados, sem grade.
* Fonte Courier Prime em todo texto da figura. O notebook carrega a fonte (licença OFL) e para com aviso claro se ela não estiver disponível.
* Exportar EPS (obrigatório para LaTeX) e PDF. Fontes incorporadas.
* Largura 80 mm (uma coluna) ou 165 mm (duas colunas). Linhas com pelo menos 0,3 pt depois da redução.
* Classes e variantes distinguidas por estilo de linha, marcador e hachura, nunca por cinza.
* Sem número de figura, legenda ou nome de arquivo dentro do arquivo. Eixos com título e rótulos legíveis. Ponto decimal e sem separador de milhar.
* Nomes figure1.eps, figure2.eps e figure3.eps.

Conteúdo:

| Figura | Largura proposta | Conteúdo |
|---|---|---|
| 1 | 80 mm | Distribuição cumulativa complementar de energia por classe, com ajuste e marca de xmin |
| 2 | 80 mm | α por classe para S0 e as variantes de sensibilidade, com a linha em 3/2 e os valores do Yang |
| 3 | 165 mm | Distribuição nula de Λ com o valor observado, e distribuições bootstrap de Δα por par |

Painéis: a legenda identifica as partes por a), b), c) em romano (epletters p2949) e figuras compostas são aceitas em coluna ou em largura total (EPLgraphics p1). As normas não tratam de letras dentro do arquivo. Decisão: letras a), b) em romano dentro das figuras compostas, iguais às da legenda. Courier Prime não é uma das fontes padrão (Times, Helvetica, Courier, Symbol, EPLgraphics p2), então fica incorporada ao EPS.

## 10. Reprodutibilidade e execução

* Sementes fixas, registradas no cabeçalho de analysis_results.md.
* Modo completo (B_BOOT, B_NULL e N_GOF nos valores padrão) gera os números do manuscrito. Modo rápido (valores reduzidos) serve para testar o notebook e grava no cabeçalho que não vale para publicação.
* O notebook estima o tempo de execução antes das etapas pesadas e permite paralelizar réplicas.
* Se o VizieR estiver inacessível, o notebook para com mensagem clara. Não há cópia local de dados.
* Idioma: analysis_results.md, os textos do notebook, os rótulos das figuras e os comentários estão em inglês. Esta especificação permanece em português. O estilo vale nos dois idiomas: sem travessões e sem ponto e vírgula. Em inglês, ortografia britânica consistente (norma da EPL, epletters p2949).

## 11. Parâmetros a aprovar

| Parâmetro | Valor proposto | Origem |
|---|---|---|
| LOGG_MIN | 4,0 (variante 4,2) | shibayama2013 p6 e p7, walkowicz2011 p4 |
| RAD_MAX | RAD_FATOR 1,5 vezes a mediana por intervalo de 200 K (variantes 1,25 e 2,0) | escolha de projeto, seção 3.2 |
| Modelo principal | M2 (truncado) | escolha de projeto a partir de D-D |
| Bootstrap | por estrela | escolha de projeto |
| B_BOOT, B_NULL, N_GOF | 1000, 500, 200 | escolha de projeto |
| Grade de xmin | 0,02 dex, com checagem | escolha de projeto |
| Nível de significância | 0,05 global, Holm nos pares | escolha de projeto |
| DELTA_EQ | 0,10 | escolha de projeto |
| N_MIN (S6) | 50 | escolha de projeto |
| Idioma de analysis_results.md e dos textos do notebook | inglês | decisão do usuário |

As demais linhas desta tabela foram aprovadas pelo usuário. RAD_FATOR e a definição por mediana de 200 K entram como parte da aprovação de RAD_MAX.

## 12. Limitações a declarar e itens sem cobertura na coleção

Limitações que o notebook e o manuscrito devem declarar:

* Sem injeção e recuperação no catálogo e sem completude por estrela (yang2019 p5).
* Energia com erro de cerca de 60 por cento e subestimada em cerca de 25 por cento pela cadência longa (yang2019 p5, yang2018 p1).
* A taxa e o expoente podem depender da rotação, com resultados de sentidos opostos entre wu2015 p1 (γ maior para rotação mais rápida) e okamoto2021 p16 a p19, Tabela 4 (índice menos inclinado para Prot menor). As classes M, K e G têm distribuições de Prot diferentes (mcquillan2014 p5, Tabela 3), de modo que uma diferença de α entre classes pode refletir mistura de rotação. A análise por rotação fica fora do P1.
* A lei de potência observada não prova criticalidade (bonachela2009 p1, watkins2016 p27).

Contagens reais da configuração principal (E7a): 3302 de 3420 estrelas casadas, 2894 estrelas e 138414 flares após o filtro de anãs, e amostra final M 549 estrelas e 44151 flares, K 1232 e 62553, G 616 e 22508, F 475 e 8938. Os números do Yang por tipo (seção 6.4) incluem gigantes e as estrelas excluídas aqui.

Itens de método sem referência na coleção, tratados como procedimentos descritos no próprio texto e não citados: bootstrap por estrela, correção de Holm, teste Q de Cochran, calibração da razão de verossimilhanças por simulação sob H0. Citar qualquer um exigiria baixar a fonte (regra R2).

## 13. Extensão com estrelas F

Motivo da extensão: o Yang conclui que as FFDs de F a M seguem lei de potência com α em torno de 2 (yang2019 p1, p8 e p19), então F pertence à pergunta de universalidade. Motivo de ficar fora do teste principal: a contaminação por pulsadores δ Scuti e γ Doradus é a principal fonte de falsos positivos do catálogo (yang2019 p3, que exclui mais de 3000 pulsadores de cerca de 8300 candidatos, p4), e na análise de flares de Van Doorsselaere a duração e a energia de F se agrupam com as das gigantes, não com G, K e M (vandoorsselaere2017 p1). Envelope convectivo existe para Teff abaixo de 6500 K (mcquillan2014 p2).

Definição da classe F: 5900 K <= Teff < 7200 K. O limite superior é F0V na escala de Pecaut (pecaut2013 p9, Tabela 5). O Yang usa 6000 a 7500 K e Van Doorsselaere 6100 a 7600 K (yang2019 p7, vandoorsselaere2017 p3).

Análises da extensão, todas rotuladas como secundárias em analysis_results.md:

1. Ajuste da classe F com o mesmo procedimento das seções 5 e 6.
2. Teste global de igualdade com quatro classes (H0F), com a calibração por simulação da seção 6.2.
3. Teste por par F contra G (H0FG), com correção de Holm junto com os outros pares.
4. Variante S10: F restrita a 5900 K <= Teff < 6500 K (F5V em 6510 K, pecaut2013 p9), que reduz a parte mais quente, onde pulsadores dominam.

O veredito principal (seção 8.2) vem do teste com M, K e G e não depende da extensão. Se a extensão discordar do teste principal, o relato descreve as duas conclusões e a lista de causas possíveis, entre elas a contaminação por pulsadores, sem escolher uma. O filtro de anãs, o RAD_MAX por intervalo de 200 K e as sensibilidades S1 a S9 valem também para F. Estrelas A (Teff acima de 7200 K) continuam fora, porque o próprio Yang encontra α diferente e sugere outro mecanismo (yang2019 p1, p7 e p8).

Saídas adicionais: linha F na tabela principal, AlphaF e demais macros com sufixo F, e a classe F na Figura 2 com estilo de linha próprio.

## 14. Análise principal com energia escalada pela luminosidade (CR-01)

Origem: Change Request CR-01, aprovado pelo usuário em 2026-10-18 (DEC-34). Ele altera as decisões D-C e D-E do DEC-24 para a análise principal e mantém a análise em energia como publicada como linha de base.

### 14.1 Motivo

O pico da distribuição diferencial de energia depende da classe em E e se alinha quando cada energia é dividida pela luminosidade da estrela. Com o pico em log E de 32,7 (M), 33,5 (K), 34,1 (G) e 34,3 (F), a dispersão é de 1,6 dex, e em log E_L é de 0,2 dex. O limiar de detecção é uma amplitude relativa (shibayama2013 p8, gunther2020 p4) e a energia é a duração equivalente vezes a luminosidade na banda Kepler (davenport2016 p5), de modo que o limite de completude em energia cresce com a luminosidade. A diferença de α entre classes em energia pura mistura completude com física. Os subsídios numéricos estão no caderno.

### 14.2 Variável

E_L = E / (L / L_sol), com L / L_sol = Rad^2 (Teff / 5772 K)^4, usando Rad e Teff do Mathur. Limites: L_bol não é a luminosidade na banda Kepler e a razão retém um fator que depende de Teff pela resposta da banda. O raio tem incerteza típica de 27 por cento (mathur2017 p1), que entra como ruído em E_L e é testada como sensibilidade (S4).

### 14.2b Parâmetros

| Parâmetro | Valor |
|---|---|
| Limite inferior comum de log E_L (XMIN_L) | 34,5 |
| Janela do ajuste truncado (WINDOW_L) | 34,5 a 35,5 |
| Janelas de sensibilidade | 34,6 a 35,4, 34,8 a 35,6, 35,0 a 36,0, 34,3 a 35,3, 34,7 a 35,7 |
| Réplicas de bootstrap por estrela, análise principal | 1000 |
| Réplicas nas variantes | 200 |
| Amostras sintéticas para a aderência com limite fixo | 200 |

### 14.3 Estimadores

* M3, lei de potência com corte exponencial acima de XMIN_L, com índice α_c e escala de corte λ_c por classe. É o modelo principal.
* M2 com limites fixos na janela WINDOW_L, estimativa secundária.
* M1 acima de XMIN_L, referência.

Todos dependem só de três somas por estrela (número de flares, soma de ln(E_L / xmin) e soma de E_L / xmin). O bootstrap reamostra estrelas inteiras, inclusive as sem flares acima do limite, e usa as somas por estrela.

### 14.4 Teste de igualdade

* H0: o índice de corte é igual em M, K e G, cada classe com sua própria escala de corte. Extensão com F como secundária.
* Estatística: soma ponderada dos desvios quadráticos dos índices em relação à média ponderada, com pesos iguais ao inverso das variâncias do bootstrap. A distribuição nula vem do bootstrap por estrela com os índices observados subtraídos, e p = (1 + número de réplicas com estatística maior ou igual à observada) / (B + 1). O p do qui-quadrado com a mesma estatística é reportado como referência.
* Diferenças por pares com intervalos de percentis e correção de Holm. Equivalência dentro da margem DELTA_EQ de 0,10 pelo intervalo de 90 por cento.
* A razão de verossimilhanças aninhada da seção 6.1 é reportada para o modelo de corte, mas trata os flares como independentes e é anticonservadora. O efeito de desenho (razão entre o erro de bootstrap por estrela e o erro de Hill) mostra quanto. Por isso o p do bootstrap por estrela é o p principal, e a calibração por simulação da seção 6.2 só vale para dados sem agrupamento.
* Aderência do modelo de corte com limite fixo: amostras sintéticas do modelo ajustado com o número observado de flares, ajuste de novo dos dois parâmetros e fração de distâncias KS sintéticas iguais ou maiores que a observada. Também trata os flares como independentes.

### 14.5 Veredito

* H0 rejeitada se o p do bootstrap por estrela for menor que 0,05. Se for rejeitada e todas as diferenças por pares forem menores que 0,10, a diferença é descrita como detectável e pequena.
* O veredito é robusto se for o mesmo em todas as variantes da seção 7 em que as três classes mantêm pelo menos 1000 flares acima do limite inferior. As variantes que mudam o veredito são listadas.
* Comparação com 3/2: rejeitada para uma classe quando 3/2 fica fora do intervalo de 99 por cento do bootstrap. Sem afirmação sobre classe de universalidade.

### 14.6 Adaptações das sensibilidades

* S1 passa a mostrar os limites inferiores fixos em energia como publicada (log E de 33,0, 33,5 e 34,0), sem escalonamento, para evidenciar o efeito da completude.
* S4 inclui, além dos fatores de 1,25, ruído de 0,2 dex na energia e ruído na luminosidade vindo dos intervalos de confiança de Rad e Teff do Mathur, sorteado uma vez por estrela em cada réplica.
* S5 passa a significar as janelas e limites de E_L listados acima. S9 é substituída por elas.
* S6 usa o índice truncado por estrela na janela de E_L, com pelo menos 50 flares na janela, o Q de Cochran e a dispersão entre estrelas de DerSimonian e Laird. As variâncias vêm da informação observada e supõem flares independentes dentro da estrela.
* S8 usa classes de Yang (4000, 5000, 6000 e 7500 K), sem filtro de anãs, bins de 0,2, 0,25 e 0,3 dex e as faixas de energia da Fig. 3 do Yang.

### 14.7 Validações acrescentadas e critério revisado

* Recuperação do índice do modelo de corte ajustado pelas somas por estrela, em amostras sintéticas.
* Tamanho do teste de igualdade com estrelas de índices diferentes em torno do valor da classe, esperado entre 0,02 e 0,10 em 0,05, e poder esperado de pelo menos 0,8 para diferença de 0,2 no índice. O mesmo teste com estatística de Wald que trata flares como independentes serve de contraste.
* Critério revisado das validações contra a lognormal (DEC-34): poder acima do tamanho do teste e crescente com o tamanho da amostra. O poder limitado do teste de aderência com xmin pelo KS é uma limitação documentada.

## 15. Etapa E8: índice por faixas estreitas de Teff

Escopo definido pelo usuário em 2026-10-19. Motivo: a igualdade entre as classes M, K e G não foi robusta na E7c, e os limites de classe do Yang e a janela de 34,3 a 35,3 invertiam o veredito. Uma explicação possível é que α varia com Teff em escala mais fina que as classes. O diagnóstico exploratório está no caderno.

### 15.1 Perguntas e hipóteses

* H0T: o índice de corte α_c é igual em todas as faixas de Teff. H1T: varia com Teff.
* Tendência: inclinação de α_c em Teff, por 1000 K, diferente de zero.
* Forma: um degrau por classe (M, K, G, F) descreve as faixas tão bem quanto as próprias faixas?

### 15.2 Faixas

* Faixas principais: 12 faixas com o mesmo número de estrelas, definidas pelos quantis de Teff das estrelas da amostra principal (M, K, G e F, 2872 estrelas), com cerca de 239 estrelas cada. Os limites das faixas são fixos nas réplicas do bootstrap.
* Variantes: 8 e 16 faixas por quantis, e faixas fixas de 250 K com fusão das faixas de cauda pequena.
* Uma faixa com menos de 1000 flares acima do limite inferior é mostrada, mas sinalizada e fica fora dos testes globais e de tendência (clauset2009 p13). No diagnóstico, a faixa mais quente tem 857 flares.

### 15.3 Estimadores

Os mesmos da seção 14.3 em cada faixa, com limite inferior comum log E_L de 34,5 e janela de 34,5 a 35,5: M3 (principal), janela truncada e M1 (referência). Bootstrap por estrela dentro de cada faixa com B de 1000 (200 nas variantes). Variante de limite relativo ao pico: limite inferior igual ao pico da faixa em log E_L (histograma de 0,2 dex) mais 0,3 dex e mais 0,5 dex, e janelas de 1 dex a partir dele. Diagnóstico obrigatório: o pico por faixa e sua inclinação em Teff, que mostra se a completude em E_L depende de Teff (no diagnóstico o pico fica entre 34,08 e 34,42, sem tendência).

### 15.4 Testes

* Igualdade em todas as faixas válidas: mesma estatística da seção 14.4, bootstrap por estrela recentrado, p = (1 + número de réplicas com estatística maior ou igual) / (B + 1).
* Tendência: regressão linear ponderada de α em Teff mediano, com pesos do inverso da variância do bootstrap, e inclinação por 1000 K com erro e intervalo vindos do bootstrap por estrela (a regressão é refeita em cada réplica).
* Forma: χ² ponderado de três modelos para α em função da faixa (constante, linear, um valor por classe), com graus de liberdade, para comparar o degrau por classe com a tendência.
* Os mesmos testes para a janela truncada e a lei de potência pura como estimativas secundárias.

### 15.5 Sistemáticas e sensibilidades

* Ruído em Teff: o Teff de cada estrela é sorteado dentro dos intervalos de confiança do Mathur em cada réplica e as estrelas são reatribuídas às faixas, o que mostra o borramento das faixas.
* Limite inferior: 34,6, 34,8 e relativo ao pico (seção 15.3).
* Exclusões: sem filtro de anãs, sem estrelas de vizinha (NS igual a p), sem estrelas sem rotação, sem estrelas não DSEP.
* A dependência de α com a rotação fica fora do P1 (BACKLOG, DEC-06). O diagnóstico mostra a mediana de Prot por faixa só como contexto, e a limitação é declarada.
* Verdade esperada de limitações: o corte exponencial tem degenerescência entre índice e escala de corte, que eleva o erro de α_c em algumas faixas (0,15 a 0,45 no diagnóstico), e a janela truncada é mais precisa, mas não é o modelo principal.

### 15.6 Veredito

H0T é rejeitada se o p do bootstrap for menor que 0,05. O veredito é robusto se for o mesmo em todas as variantes com pelo menos 1000 flares na cauda das faixas usadas. A descrição da forma de α(Teff) é qualitativa e só se afirma o que a tabela e o teste de forma sustentarem.

### 15.7 Saídas

* Novas seções em analysis_results.md: faixas, α por faixa para os três estimadores, diagnóstico do pico, testes, formas, variantes e limitações automáticas.
* Novos macros em results_macros.tex: número de faixas, p de igualdade, inclinação e erro por estimador, p da forma, pico e sua inclinação.
* Figura 4 (80 mm, mesmo estilo e fonte): α_c e α da janela contra Teff mediano, com erro de bootstrap, linhas verticais nos limites de classe (3850, 5300 e 5900 K), linha horizontal em 3/2 e símbolos abertos para faixas com cauda abaixo de 1000. O limite de 7 páginas da EPL exige escolher no manuscrito quais figuras ficam, o que não é decidido aqui.

### 15.8 Validações

Em dados sintéticos com estrelas heterogêneas e α dependente da faixa: o teste global mantém o tamanho (0,02 a 0,10 em 0,05) sem tendência, detecta uma tendência de 0,1 por 1000 K com poder de pelo menos 0,8, e a inclinação é recuperada sem viés maior que 0,02 por 1000 K.

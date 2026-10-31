# Calculadora de Propagação de Incertezas — Frequência de um Fio Esticado

Calcula a **frequência de oscilação** de um fio esticado (ondas estacionárias) e a sua **incerteza**, propagando os erros das medidas feitas no laboratório. Tem **interface gráfica**, modo **terminal** e uma versão **web** (arquivo HTML único).

## Modelo físico

$$
f = \frac{n}{2L}\sqrt{\frac{T}{\mu}}, \qquad T = m_{obj}\,g, \qquad \mu = \frac{m_{fio}}{l_{fio}}
$$

| Símbolo | Significado | Unidade |
|---|---|---|
| $f$ | frequência | Hz |
| $n$ | número de ventres (exato, sem incerteza) | — |
| $L$ | comprimento vibrante do fio | m |
| $m_{obj}$ | massa suspensa que gera a tensão | kg |
| $m_{fio}$ | massa do fio | kg |
| $l_{fio}$ | comprimento total do fio (usado para $\mu$) | m |
| $g$ | aceleração da gravidade local | m/s² |

### Propagação de incertezas

Com as grandezas independentes, a incerteza de $f$ vem da propagação linear (derivadas parciais):

$$
\sigma_f^2 = \left(\frac{\partial f}{\partial L}\sigma_L\right)^2 + \left(\frac{\partial f}{\partial T}\sigma_T\right)^2 + \left(\frac{\partial f}{\partial \mu}\sigma_\mu\right)^2
$$

que equivale, em termos relativos, a:

$$
\left(\frac{\sigma_f}{f}\right)^2 = \left(\frac{\sigma_L}{L}\right)^2 + \left(\frac{\sigma_T}{2T}\right)^2 + \left(\frac{\sigma_\mu}{2\mu}\right)^2
$$

com

$$
\sigma_T = T\sqrt{\left(\frac{\sigma_{m_{obj}}}{m_{obj}}\right)^2 + \left(\frac{\sigma_g}{g}\right)^2}, \qquad
\sigma_\mu = \mu\sqrt{\left(\frac{\sigma_{m_{fio}}}{m_{fio}}\right)^2 + \left(\frac{\sigma_{l_{fio}}}{l_{fio}}\right)^2}
$$

### Hipóteses

- Fio ideal, sem atrito na polia: a tensão é o peso do objeto suspenso.
- $n$ é exato.
- Grandezas medidas **independentes** (sem correlação entre elas).
- Incertezas pequenas o bastante para a aproximação linear valer (o Monte Carlo ajuda a conferir).

## Funcionalidades

- Interface gráfica (Tkinter) que atualiza o resultado enquanto você digita.
- Aceita vírgula ou ponto decimal.
- Escolha de unidades: comprimentos em m, cm, mm e massas em kg, g, mg (conversão automática para SI).
- Resultado com arredondamento correto: σ com 2 algarismos significativos e o valor no mesmo decimal.
- Contribuição de cada variável (L, T, μ) para a incerteza, para saber o que medir melhor.
- Verificação por **Monte Carlo** independente do cálculo analítico.
- Comparação com a frequência medida pelo **erro normalizado** $E_n = \dfrac{|f - f_{exp}|}{\sqrt{\sigma_f^2 + \sigma_{f_{exp}}^2}}$ (compatível se $|E_n| \le 1$).
- Validação das entradas (valores > 0, incertezas ≥ 0).

## Requisitos

- Python 3.9 ou superior.
- **Nenhuma biblioteca externa** (usa só a biblioteca padrão).
- Para a interface gráfica, o Tkinter (já vem com o Python no Windows e no macOS). No Ubuntu/Debian:

```bash
sudo apt install python3-tk
```

Se o Tkinter não estiver disponível, o programa cai automaticamente para o modo terminal.

## Como usar

### Interface gráfica (padrão)

```bash
python propagacao_incertezas.py
```

Preencha valores e incertezas, escolha as unidades e veja o resultado. Os botões **Exemplo**, **Monte Carlo** e **Copiar resultado** ficam abaixo dos campos.

### Terminal interativo

```bash
python propagacao_incertezas.py --cli
```

O programa pergunta cada valor.

### Linha de comando com argumentos

```bash
python propagacao_incertezas.py \
  --n 3 --L 80 --sL 0.2 \
  --m-obj 200 --sm-obj 1 \
  --m-fio 1.2 --sm-fio 0.05 \
  --l-fio 200 --sl-fio 0.5 \
  --g 9.78 --mc --f-exp 107 --sf-exp 1
```

| Argumento | Descrição |
|---|---|
| `--n` | número de ventres |
| `--L`, `--sL` | comprimento vibrante e sua incerteza |
| `--m-obj`, `--sm-obj` | massa suspensa e sua incerteza |
| `--m-fio`, `--sm-fio` | massa do fio e sua incerteza |
| `--l-fio`, `--sl-fio` | comprimento total do fio e sua incerteza |
| `--g`, `--sg` | gravidade local (padrão 9,78 m/s²) e sua incerteza (padrão 0) |
| `--uc` | unidade dos comprimentos: `m`, `cm` (padrão) ou `mm` |
| `--um` | unidade das massas: `kg`, `g` (padrão) ou `mg` |
| `--f-exp`, `--sf-exp` | frequência medida (Hz) e sua incerteza, para comparação |
| `--mc` | roda a verificação por Monte Carlo |
| `--cli` | força o modo terminal |

### Saída de exemplo

```
Frequência f = (107.1 ± 2.3) Hz   [2.12 %]
Tensão T = 1.9560 ± 0.0098 N
Densidade linear μ = 0.6000 ± 0.0250 g/m

Contribuição para a variância de f:
  mu (m_fio e l_fio)       97.2 %
  L                         1.4 %
  T (m_obj e g)             1.4 %

Monte Carlo: f = 107.1278 Hz, σ_f = 2.2771 Hz (analítico: 2.2662 Hz)

Comparação com f_exp = 107.0 Hz: E_n = 0.02 → compatível (|E_n| ≤ 1)
```

Neste exemplo a massa do fio domina a incerteza (97 %): medir um fio mais longo ou usar uma balança de melhor resolução reduz bastante o $\sigma_f$.

### Como biblioteca

```python
from propagacao_incertezas import Entrada, calcular

e = Entrada(n=3, L=0.80, sL=0.002, m_obj=0.200, sm_obj=0.001,
            m_fio=0.0012, sm_fio=0.00005, l_fio=2.0, sl_fio=0.005, g=9.78)
r = calcular(e)
print(r.f, r.sigma_f)
```

Todas as grandezas de `Entrada` são em unidades SI (m, kg, m/s²).

## Versão web

O arquivo `calculadora-incertezas.html` é a mesma calculadora em uma página única, que abre direto no navegador, sem internet e sem Python. Com o GitHub Pages ativado no repositório, ele vira um link público.

## Estrutura

```
.
├── propagacao_incertezas.py     # cálculo + interface gráfica + terminal
├── calculadora-incertezas.html  # versão web (arquivo único)
└── README.md
```

## Licença

Projeto livre para uso acadêmico e educacional.

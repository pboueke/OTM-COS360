# OTM-COS360

## Trabalho da disciplina [Otimização (COS360) - UFRJ](http://www.cos.ufrj.br/~luidi/cos360/otim.html)

**2016/P2**: Terças e quintas de 15:00-17:00

**Professor**: Luidi Simonetti

**Alunos**: Guilherme Thurler e Pedro Boueke

**[Tarefa](https://github.com/pboueke/OTM-COS360/blob/master/doc/COS360_Trabalho_09_2016.pdf)**

---

### Funções

1. f(x1, x2) = ln(1 + ln(x1)^2 + ln(x2)^2)
2. f(x1, x2) = ln(1 + x1^2 + (x1^2 - x2)^2)

---

#### Tarefa 1 - *Faça um pequeno estudo para um melhor entendimento do comportamento de suas funções.*

#### f(x1, x2) = ln(1 + ln(x1)^2 + ln(x2)^2)

* **[Gráfico](https://www.google.com.br/search?safe=off&q=ln%281+%2B+ln%28x%29%5E2+%2B+ln%28y%29%5E2%29&oq=ln%281+%2B+ln%28x%29%5E2+%2B+ln%28y%29%5E2%29&gs_l=serp.3..0i8i30k1.3682.6299.0.7751.4.4.0.0.0.0.305.861.0j3j0j1.4.0....0...1c.1.64.serp..0.3.666.kWqgk3Mgvus)**
![function 1](https://github.com/pboueke/OTM-COS360/blob/master/img/func1.gif?raw=true | width=250)

* **Derivada**
##### (2\*ln(x1)) / x1\*(1 + ln(x1)^2 + ln(x2)^2) em x1
##### (2\*ln(x2)) / x2\*(1 + ln(x2)^2 + ln(x1)^2) em x2

* **Mínimo**: (1, 1)

#### f(x1, x2) = ln(1 + x1^2 + (x1^2 - x2)^2)

* **[Gráfico](https://www.google.com.br/search?safe=off&q=ln%281+%2B+x%5E2+%2B+%28x%5E2+-+y%29%5E2%29&oq=ln%281+%2B+x%5E2+%2B+%28x%5E2+-+y%29%5E2%29&gs_l=serp.3..0i8i30k1l10.32603.39538.0.39688.6.6.0.0.0.0.254.790.0j3j1.4.0....0...1c.1.64.serp..2.2.364...0i8i10i30k1.8JgrcN1a4-M)**
![function 2](https://github.com/pboueke/OTM-COS360/blob/master/img/func2.gif?raw=true | width=250)

* **Derivada**
##### (4\*x1\*(x1^2 + x2) + 2\*x1) / ((x1^2 + x2)^2 + x1^2 + 1) em x1
##### (2\*(x2 + x1^2)) / ((x2 + x1^2)^2 + x1^2 + 1) em x2

* **Mínimo**: (0, 0)

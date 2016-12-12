# OTM-COS360

## Relatório do trabalho da disciplina [Otimização (COS360) - UFRJ](http://www.cos.ufrj.br/~luidi/cos360/otim.html)

**2016/P2**: Terças e quintas de 15:00-17:00

**Professor**: Luidi Simonetti

**Alunos**: Guilherme Thurler e Pedro Boueke

**[Tarefa](https://github.com/pboueke/OTM-COS360/blob/master/doc/COS360_Trabalho_09_2016.pdf)**

---

### Funções

1. f(x1, x2) = ln(1 + ln(x1)^2 + ln(x2)^2)
2. f(x1, x2) = ln(1 + x1^2 + (x1^2 + x2)^2)

---

#### Tarefa 1 - *Faça um pequeno estudo para um melhor entendimento do comportamento de suas funções.*

#### f(x1, x2) = ln(1 + ln(x1)^2 + ln(x2)^2)

* **Gráfico**
![function 1](https://github.com/pboueke/OTM-COS360/blob/master/img/func1.gif?raw=true | width=250)

* **Derivada**
##### (2\*ln(x1)) / x1\*(1 + ln(x1)^2 + ln(x2)^2)

* **Mínimo**: (1, 1)

#### f(x1, x2) = ln(1 + x1^2 + (x1^2 + x2)^2)

* **Gráfico**
![function 2](https://github.com/pboueke/OTM-COS360/blob/master/img/func2.gif?raw=true | width=250)

* **Derivada**
##### (4\*x1\*(x1^2 + x2) + 2\*x1) / ((x1^2 + x2) + x1^2 + 1)

* **Mínimo**: (0, 0)

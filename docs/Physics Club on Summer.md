---
tags: [maths, physics, physics club, school]
title: Physics Club on Summer
created: '2022-06-22T13:52:45.472Z'
modified: '2022-06-22T16:44:54.665Z'
---

# Physics Club on Summer
## Preliminari di Matematica

Qualche strumento matematico utile per i problemi di Fisica

### Derivate

#### Approssimazioni
Per una funzione $f(x)$ cerchiamo una approssimazione per $x\sim x_0$

$$
x = x_0 + \Delta x \qquad \textrm{ con } \qquad \Delta x \sim 0
$$ 

L'approssimazione della forma più semplice ha l'aspetto di un polinomio:

$$
f(x) = f(x_0+\Delta x) \approx a_0 + a_1 (x-x_0) + a_2 (x-x_0)^2 + \dots\\
\ \\
= a_0 + a_1 \Delta x + a_2 \Delta x^2 + \dots
$$

I puntini di sospensione rappresentano dei termini ulteriori che non esplicitiamo, se tutto va bene possiamo trascurarli se $\Delta x$ è _"abbastanza piccolo"_ (termine da specificare). In generale al posto dei puntini di sospensione avremmo un termine con una potenza superiore di $\Delta x$, per esempio nel caso di una approssimazione di ordine 2 (con un polinomio di secondo grado) abbiamo:

$$
f(x) = f(x_0+\Delta x) \approx a_0 + a_1 (x-x_0) + a_2 (x-x_0)^2 + R_3(x_0,x)
$$

laddove il termine $R_3(x_0,x)$ deve rispettare queste condizioni:

$$ 
x\rightarrow x_0 \Rightarrow\\
\ \\
\qquad\qquad\qquad\Rightarrow R_3(x_0,x) \rightarrow 0\\
\qquad\qquad\qquad\Rightarrow \frac{R_3(x_0,x)}{x-x_0} \rightarrow 0\\
\qquad\qquad\qquad\Rightarrow \frac{R_3(x_0,x)}{(x-x_0)^2} \rightarrow 0\\
$$
#### Approssimazioni di grado 0
Una funzione $f(x)$ si dice **continua** in $x=x_0$ se possiamo scrivere l'equazione seguente, che ci dà una approssimazione di $f(x)$ con una costante (vedi lezione 1)

$$
x\sim x_0 \qquad f(x) = f(x_0) + R_1(x_0, x)\\
\ \\
\qquad\qquad\qquad\qquad\textrm{ dove } \qquad x\rightarrow x_0 \Rightarrow R_1(x_0, x) \rightarrow 0
$$

#### Approssimazioni lineari (di grado 1)
Una funzione $f(x)$ si dice **derivabile** in $x=x_0$ se possiamo scrivere l'equazione seguente, che ci dà una approssimazione di $f(x)$ con una retta (vedi lezione 1)

$$
x\sim x_0 \qquad f(x) = f(x_0) + m\cdot(x-x_0) + R_2(x_0, x)\\
\ \\
\qquad\qquad\qquad\qquad\textrm{ dove } \qquad x\rightarrow x_0 \Rightarrow R_2(x_0, x) \rightarrow 0\\
\ \\
\qquad\qquad\qquad\qquad\phantom{\textrm{ dove }} \qquad \qquad \qquad \Rightarrow \frac{R_2(x_0, x)}{x-x_0} \rightarrow 0
$$

In questa equazione la pendenza $m$ si chiama **derivata di** $f(x)$ in $x=x_0$ cioè

$$
m = f'(x_0)
$$

in notazione usuale. Ricaviamo anche la forma alternativa

$$
\frac{f(x_0+\Delta x) - f(x_0)}{\Delta x} \rightarrow f'(x_0) \qquad \textrm{ per } \Delta x \rightarrow 0
$$

##### Interpretazione geometrica (retta tangente)
La nostra approssimazione lineare ha una diretta interpretazione geometrica: se prendiamo i termini fino al primo ordine in $\Delta x$:

$$
x\sim x_0 \qquad f(x) \approx f(x_0) + m\cdot(x-x_0)\\
\ \\ 
\qquad\qquad\qquad\qquad y = f(x_0) + f'(x_0)\cdot(x-x_0)
$$

otteniamo l'equazione della retta tangente al grafico della funzione $f(x)$ nel punto $x=x_0$.

## Regole di derivazione
Calcolare la pendenza della retta tangente, ovvero il parametro della approssimazione lineare, è un compito che può venir portato a termine una volta data la funzione: in effetti si tratta di seguire un algoritmo ben determinato, algoritmo che stiamo imparando poco per volta nelle nostre lezioni.

### Linearità della derivata
Per _linearità_ si intendono queste due relazioni (le funzioni $f(x)$ e $g(x)$ sono entrambe _derivabili_ - vedi definizione qui sopra):
$$
\begin{array}{rclcrcl}
h(x) & = & f(x) + g(x) & \Rightarrow & h'(x) & = & f'(x) + g'(x)\\
h(x) & = & c \cdot f(x) & \Rightarrow & h'(x) & = & c\cdot f'(x) \qquad\qquad (c\in\mathbb{R})
\end{array}
$$
### Derivata di Polinomi
Nelle lezioni abbiamo visto che, se $k\in\mathbb Z$
$$
\begin{array}{rclcrcl}
h(x) & = & x^k & \Rightarrow & h'(x) & = & k \cdot x^{k-1}
\end{array}
$$
Insieme alle [proprietà di linearità](#linearità-della-derivata) possiamo usare questa relazione per determinare la derivata di qualunque polinomio:
$$
h(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_N x^N \Rightarrow h'(x) = a_1 + 2 a_2 x + 3 a_3 x^2 + \cdots + N a_N x^{N-1}
$$
per esempio
$$
h(x) = 3 x + 5 x^3 - 4 x^5 + x^7 \Rightarrow h'(x) = 3 + 15 x^2 - 20 x^4 + 7 x^6
$$

## Equazioni differenziali
Una **Equazione Differenziale Ordinaria** (**EDO**) è una equazione in cui l'incognita è una **funzione** della variabile indipendente (per esempio $y(x)$), e in questa equazione compare la derivata (o _le_ derivate) della funzione incognita. Per esempio nella lezione 3 e nella lezione 4 abbiamo visto delle equazioni differenziali a _integrazione immediata_. Per esempio:
$$
y'(x) = x^2
$$
Una equazione di questo genere ammette in genere una famiglia infinita di soluzioni, dipendenti da parametri che si chiamano *costanti di integrazione*. Per esempio l'equazione qui sopra ammette come soluzioni le funzioni della famiglia di equazioni
$$
\left\{ y(x) = \frac{1}{3} x^3 + c \ , \ c\in\mathbb{R} \right\}
$$
### Problemi di Cauchy
Nella pratica risolviamo problemi di equazioni differenziali in cui oltre all'equazione differenziale propriamente detta abbiamo delle altre condizioni _accessorie_, che tipicamente riguardano il valore della funzione incognita per un qualche valore iniziale della funzione, per esempio un problema di Cauchy potrebbe essere (tratto dall'EDO qui sopra)
$$
\left\{
  \begin{array}{rcl}
    y'(x) & = & x^2 \\
    y(0) & = & 3 \\
  \end{array}
\right.
$$
La condizione iniziale ci permette di determinare, all'interno della famiglia di _soluzioni generali_ $\left\{ y(x) = \frac{1}{3} x^3 + c \ , \ c\in\mathbb{R} \right\}$ la soluzione particolare di questo problema di Cauchy, tramite sostituzione:
$$
3 = y(0) = \frac{1}{3} 0^3 + c \Rightarrow c = 3
$$
Cioè soluzione del problema di Cauchy sarà la funzione
$$
y(x) = \frac{1}{3}x^3 + 3
$$
che possiamo verificare rispetta sia l'equazione differenziale $y'(x) = x^2$ sia la condizione iniziale $y(0)=3$.
$$
y'(x) = \frac{1}{3}(x^3)' + (3)' = \frac{1}{3} 3 x^2 + 0 = x^2 \qquad y(0) = \frac{1}{3} 0^3 + 3 = 3
$$



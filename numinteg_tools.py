import numpy as np
import matplotlib.pyplot as plt

def euler(x0, ts, f):
    """Algoritmo di Euler in avanti.
    
        euler(x0, ts, f) -> restituisce un array numpy, con tante colonne quante sono le componenti 
                            dell'array x(t), e una riga per ogni tempo t_k, con i valori previsti per
                            le quantità x(t) regolate dall'equazione differenziale
                                    dx(t) / dt = f(x, t)
        
        x0 è la condizione iniziale (array o lista di valori iniziali); il numero di componenti di x0
           viene usato per dedurre il numero di equazioni del sistema.
        ts array numpy contenente i tempi della discretizzazione; la lunghezza di questo array sarà
           uguale alla lunghezza della simulazione.
        f  funzione che dà la velocità di variazione delle quantità x(t); deve essere una funzione che 
           accetta due variabili, x e t, la prima un array o lista di valori delle componenti di x(t) 
           (quindi con la stessa lunghezza del vettore delle condizioni iniziali), e la seconda uguale 
           al tempo al quale viene calcolata la velocità.
           """
    
    # inizializza l'array dei valori calcolati
    xs = np.zeros(shape=(len(ts), len(x0)), dtype=np.float)
    
    # condizione iniziale -> primo elemento dell'array
    xs[0] = x0     
    
    i = 0
    while i<len(ts)-1:         # ciclo su tutti i tempi
        dt = ts[i+1] - ts[i]   # ampiezza dell'intervallo di tempo (in linea di principio può essere non costante)
        v = f(xs[i], ts[i])    # velocità di variazione: è un vettore calcolato 
        
        xs[i+1] = xs[i] + dt * v    # il passo di update che permette di calcolare il valore successivo
                                    # usando la velocità di variazione calcolata sopra
        
        i += 1
        
    return xs
    
   
def rk2(x0, ts, f, α=1, β=1, a=0.5):
    """Algoritmo di Runge-Kutta a 2 stadi (esplicito).
    
        rk2(x0, ts, f) -> restituisce un array numpy, con tante colonne quante sono le componenti 
                          dell'array x(t), e una riga per ogni tempo t_k, con i valori previsti per
                          le quantità x(t) regolate dall'equazione differenziale
                                  dx(t) / dt = f(x, t)
        
        x0 è la condizione iniziale (array o lista di valori iniziali); il numero di componenti di x0
           viene usato per dedurre il numero di equazioni del sistema.
        ts array numpy contenente i tempi della discretizzazione; la lunghezza di questo array sarà
           uguale alla lunghezza della simulazione.
        f  funzione che dà la velocità di variazione delle quantità x(t); deve essere una funzione che 
           accetta due variabili, x e t, la prima un array o lista di valori delle componenti di x(t) 
           (quindi con la stessa lunghezza del vettore delle condizioni iniziali), e la seconda uguale 
           al tempo al quale viene calcolata la velocità.
           """
    
    # inizializza l'array dei valori calcolati
    xs = np.zeros(shape=(len(ts), len(x0)), dtype=np.float)
    
    # condizione iniziale -> primo elemento dell'array
    xs[0] = x0
    
    i = 0
    while i<len(ts)-1:              # ciclo su tutti i tempi
        dt = ts[i+1] - ts[i]        # ampiezza dell'intervallo di tempo (in linea di principio può essere non costante)
        
        # primo stadio (coincide con Euler)
        v1 = f(xs[i], ts[i])      
        
        # secondo stadio: usiamo un punto intermedio
        v2 = f(xs[i] + α * dt * v1, ts[i] + β * dt)
        
        xs[i+1] = xs[i] + dt * (a * v1 + (1-a) * v2)   # il passo di update che permette di calcolare il valore successivo
                                                       # usando le velocità di variazione del primo e del secondo
                                                       # stadio calcolate sopra
        
        i += 1
        
    return xs
    

class Fields:
    @staticmethod
    def oscilloc(α=1.0, ω=2*np.pi/(24*3600)):
        """oscilloc(α, ω) genera una funzione che restituisce la velocità di variazione 
        di una popolazione di batteri secondo il modello visto sopra.
        Il parametro 
        Per avere il modello privo di oscillazioni basta porre ω=0 nell'invocare la funzione."""
        def wrap(x, t):
            return np.array([α*(1-0.5*np.sin(ω*t))*x[0]])
        return wrap

    @staticmethod
    def harmonic(ω):
        def wrap(x, t):
            return np.array([x[1], - ω**2 * x[0]])

        return wrap

    @staticmethod
    def damped_oscillator(ω, β):
        def wrap(x, t):
            return np.array([x[1], - ω**2 * x[0] - β * x[1]])

        return wrap

    @staticmethod
    def SIR(a=0.1, b=0.1):
        def g(x, t):
            s, i, r = x
            sp = -a*s*i
            ip = a*s*i - b*i
            rp = b*i

            return np.array([sp, ip, rp])

        return g

    @staticmethod
    def van_der_Pol_oscillator(ω=1.0, β=0.0):
        def g(x, t):
            return np.array([x[1], -ω**2*x[0]-β*(x[0]**2-1)*x[1]])

        return g

    @staticmethod
    def free_fall(g_x=0.0, g_y=0.0):
        def wrap(x, t):
            # x[0], x[1] x, v_x
            # x[2], x[3] y, v_y
            return np.array([x[1],g_x,x[3],g_y])

        return wrap

    @staticmethod
    def budworm(r, k):
        def wrap(x, t):
            return np.array([r*x*(1-x/k)-x**2/(1+x**2)])

        return wrap

    @staticmethod
    def bead_rotating_hoop(ϵ=0.5, γ=0.5, **kw):
        def wrap(x, t):
            return np.array([x[1], - (x[1] + np.sin(x[0]) - γ * np.sin(x[0]) * np.cos(x[0]))/ϵ])

        return wrap
    
    
def symplectic(x0, p0, ts, fx, fp):
    """Algoritmo di integrazione simplettica (I ord.).
    
        symplectic(x0, p0, ts, fx, fp) -> restituisce un array numpy, con tante colonne quante sono le componenti 
                                          dell'array x(t), e una riga per ogni tempo t_k, con i valori previsti per
                                          le quantità x(t) regolate dall'equazione differenziale
                                          dX(t) / dt = f(X, t)
                                          In questo caso X = (x, p) è formato da coordinate canoniche
                                          e impulsi coniugati
        
        x0 è la condizione iniziale (array o lista di valori iniziali) per le variabili posizione; il numero 
           di componenti di x0 viene usato per dedurre il numero di equazioni del sistema per le variabili posizione. 
        p0 è la condizione iniziale (array o lista di valori iniziali) per le variabili posizione; il numero 
           di componenti di p0 viene usato per dedurre il numero di equazioni del sistema per le variabili impulso. 
        ts array numpy contenente i tempi della discretizzazione; la lunghezza di questo array sarà
           uguale alla lunghezza della simulazione.
        fx funzione che dà la velocità di variazione delle variabili posizione x(t); deve essere una 
           funzione che accetta tre variabili, x p e t, la prima un array o lista di valori delle componenti 
           di x(t) (quindi con la stessa lunghezza del vettore delle condizioni iniziali), la seconda un array 
           o lista di valori delle componenti di p(t)e la terza uguale al tempo al quale viene calcolata la velocità.
        fp funzione che dà la velocità di variazione delle variabili quantità di moto p(t); deve essere una 
           funzione che accetta tre variabili, x p e t, la prima un array o lista di valori delle componenti 
           di x(t) (quindi con la stessa lunghezza del vettore delle condizioni iniziali), la seconda un array 
           o lista di valori delle componenti di p(t)e la terza uguale al tempo al quale viene calcolata la velocità.
           """
        
    # inizializza l'array dei valori calcolati per le x
    xs = np.zeros(shape=(len(ts), len(x0)), dtype=np.float)
    
    # inizializza l'array dei valori calcolati per le p
    ps = np.zeros(shape=(len(ts), len(p0)), dtype=np.float)
    
    # condizione iniziale -> primo elemento dell'array
    xs[0] = x0
    ps[0] = p0
    
    i = 0
    while i<len(ts)-1:                       # ciclo sui tempi
        dt = ts[i+1] - ts[i]                 # ampiezza dell'intervallo di tempo (in linea di principio può essere non costante)
        
        # primo stadio, calcola l'impulso ausiliario al tempo intermedio
        v1 = fp(xs[i], ps[i], ts[i])

        pstar = ps[i] + v1 * dt / 2   # aggiornamento
        
        # secondo stadio, determina l'update della posizione anche in funzione dell'impulso ausiliario
        v2 = fx(xs[i], pstar, ts[i])
        xs[i+1] = xs[i] + v2 * dt
        
        # ultimo stadio, determina l'update definitivo della quantità di moto
        v3 = fp(xs[i+1], pstar, ts[i])
        ps[i+1] = pstar + v3 * dt / 2 
        
        i += 1
        
    return xs, ps


class HamiltonianField:
    @staticmethod
    def double_spring(k_1, k_2, m_1, m_2):
        
        m = np.array([m_1, m_2])
        def wrap_x(x, p, t):
            return np.array(p/m)
        
        a = np.array([[-(k_1+k_2)/m_1, k_2/m_1], [k_2/m_2, -k_2/m_2 ]])
        def wrap_p(x, p, t):
            return a.dot(x)


        return wrap_x, wrap_p
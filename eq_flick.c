#include <stdio.h>
#include <stdlib.h>

// Funzione che integra l'equazione di flick usando Eulero

void flick_equation_1d(double *rho, double *D,double **rho_evol, int nx, int nt, double dx, double dt,int b) {
    // Array temporaneo per memorizzare i valori aggiornati
    double *rhon = (double *)malloc(nx * sizeof(double));
    double *Dn= (double *)malloc(nx * sizeof(double)); //utile solo se ho D che varia nel tempo.
    int n=1,c=b;
    if (Dn == NULL || rhon==NULL) {
        fprintf(stderr, "Errore di allocazione della memoria.\n");
        exit(1);
    }
    
    // Loop sui passi temporali
    for (int t = 1; t <= nt; t++) {
        //riempio array temporaneo
        for (int i = 0; i < nx; i++){
            rhon[i] = rho[i];
            Dn[i]=D[i];
        } 

        // Loop sui punti spaziali (escludendo i bordi)
                                                                            
        for (int i = 1; i < nx - 1; i++)  {
            //0.375 viene da 3/2 *1/2 *1/2 il primo viene dell'eq, gli altri due delle derivate 
            rho[i] += dt / (dx * dx) * (0.375*(Dn[i+1]-Dn[i-1])*(rhon[i+1]-rhon[i-1] )+D[i]*(rhon[i+1]+rhon[i-1]-2*rhon[i])+0.5*rhon[i]*(Dn[i+1]+Dn[i-1]-2*Dn[i]));

            if(t==c)rho_evol[n][i]=rho[i];
        }
        if(t==c){
            n++;
            c*=b;
        }
    }
    free(Dn);
    free(rhon);
}


// Funzione che integra l'equazione di flick usando Eulero in 2D
void flick_equation_2d(double *rho, double *D, int nx, int ny, int nt, double dx, double dy, double dt) {
    // Array temporaneo per memorizzare i valori aggiornati
    double *rhon = (double *)malloc(nx * ny * sizeof(double));
    double *Dn = (double *)malloc(nx * ny * sizeof(double));
    if (Dn == NULL || rhon == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria.\n");
        exit(1);
    }

    // Loop sui passi temporali
    for (int t = 0; t < nt; t++) {
        // Copia i valori attuali negli array temporanei
        for (int i = 0; i < nx * ny; i++) {
            rhon[i] = rho[i];
            Dn[i] = D[i];
        }

        // Loop sui punti spaziali (escludendo i bordi)
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                int idx = i * ny + j;             // Indice dell'elemento attuale
                int idx_ip1 = (i + 1) * ny + j;   // Indice dell'elemento a destra
                int idx_im1 = (i - 1) * ny + j;   // Indice dell'elemento a sinistra
                int idx_jp1 = i * ny + (j + 1);   // Indice dell'elemento in alto
                int idx_jm1 = i * ny + (j - 1);   // Indice dell'elemento in basso

                // Calcoli delle derivate finite e degli aggiornamenti
                double dDx = (Dn[idx_ip1] - Dn[idx_im1]) * (rhon[idx_ip1] - rhon[idx_im1]);
                double dDy = (Dn[idx_jp1] - Dn[idx_jm1]) * (rhon[idx_jp1] - rhon[idx_jm1]);
                double Dx = D[idx] * (rhon[idx_ip1] + rhon[idx_im1] - 2 * rhon[idx]);
                double Dy = D[idx] * (rhon[idx_jp1] + rhon[idx_jm1] - 2 * rhon[idx]);
                double rhoD = rhon[idx] * (Dn[idx_ip1] + Dn[idx_im1] + Dn[idx_jp1] + Dn[idx_jm1] - 4 * Dn[idx]);

                // Aggiornamento della densità rho
                rho[idx] += dt * (0.375 * (dDx + dDy) + Dx + Dy + 0.5 * rhoD) / (dx * dx);
            }
        }
    }
    // Libera la memoria degli array temporanei
    free(Dn);
    free(rhon);
}


// Funzione che integra l'equazione di flick usando Eulero in 2D con PBC
void flick_equation_2d_pbc(double *rho, double *D, int nx, int ny, int nt, double dx, double dy, double dt) {
    // Array temporaneo per memorizzare i valori aggiornati
    double *rhon = (double *)malloc(nx * ny * sizeof(double));
    double *Dn = (double *)malloc(nx * ny * sizeof(double));
    if (Dn == NULL || rhon == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria.\n");
        exit(1);
    }

    // Loop sui passi temporali
    for (int t = 0; t < nt; t++) {
        // Copia i valori attuali negli array temporanei
        for (int i = 0; i < nx * ny; i++) {
            rhon[i] = rho[i];
            Dn[i] = D[i];
        }

        // Loop sui punti spaziali
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                int idx = i * ny + j;               // Indice dell'elemento attuale
                int idx_ip1 = ((i + 1) % nx) * ny + j; // Indice dell'elemento a destra con PBC
                int idx_im1 = ((i - 1 + nx) % nx) * ny + j; // Indice dell'elemento a sinistra con PBC
                int idx_jp1 = i * ny + ((j + 1) % ny); // Indice dell'elemento in alto con PBC
                int idx_jm1 = i * ny + ((j - 1 + ny) % ny); // Indice dell'elemento in basso con PBC

                // Calcoli delle derivate finite e degli aggiornamenti
                double dDx = (Dn[idx_ip1] - Dn[idx_im1]) * (rhon[idx_ip1] - rhon[idx_im1]);
                double dDy = (Dn[idx_jp1] - Dn[idx_jm1]) * (rhon[idx_jp1] - rhon[idx_jm1]);
                double Dx = D[idx] * (rhon[idx_ip1] + rhon[idx_im1] - 2 * rhon[idx]);
                double Dy = D[idx] * (rhon[idx_jp1] + rhon[idx_jm1] - 2 * rhon[idx]);
                double rhoD = rhon[idx] * (Dn[idx_ip1] + Dn[idx_im1] + Dn[idx_jp1] + Dn[idx_jm1] - 4 * Dn[idx]);

                // Aggiornamento della densità rho
                rho[idx] += dt * (0.375 * (dDx + dDy) + Dx + Dy + 0.5 * rhoD) / (dx * dx);
            }
        }
    }
    // Libera la memoria degli array temporanei
    free(Dn);
    free(rhon);
}



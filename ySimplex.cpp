// Algorithme dual du simplexe (méthode de Lemke) (fonction simplexe_dual).

// a : le tableau des coefficients du problème ;
// sol : un tableau renvoyant la solution optimale ;
// ineq1 : le nombre des inéquations en <= ;
// ineq2 : le nombre des inéquations en >= ;
// eq : le nombre des équations ;
// n : le nombre de variables.

// Elle calcule dans le tableau sol une solution optimale du programme linéaire
// (l'opposé de la valeur correspondante de la fonction économique est calculé dans a[0][0]).

// Les données doivent être fournies dans le tableau a de la façon suivante:

//     * ligne 0 de a : la fonction économique
//       L'élément 0 est égal à 0, les éléments 1 à n aux coefficients des variables dans la fonction

//     * lignes 1 à ineq1 de a : les inéquations en <=
//       L'élément 0 est égal au second membre,
//       les éléments 1 à n aux coefficients des variables dans l'inéquation

//     * lignes ineq1+1 à ineq1+ineq2 de a : les inéquations en >=
//       L'élément 0 est égal au second membre,
//       les éléments 1 à n aux coefficients des variables dans l'inéquation

//     * lignes ineq1+ineq2+1 à ineq1+ineq2+eq de a : les équations
//       L'élément 0 est égal au second membre,
//       les éléments 1 à n aux coefficients des variables dans l'équation


#include "simplex.H"

int Simplex::ads_sortant(int m,int n,int phase)
{
    int i,j,k,l;
    double d,s,min;
    l=0;
    min=0.0;
    k=0;
    if(phase==1)
    {
        for(j=1;j<=n;j++)
        {
            if(hb[j]==n+m)
                k=j;
        }
    }
    if((phase==2)||(k!=0))
    {
        for(i=1;i<=m;i++)
        {
            d=a[i][k];
            s=0.0;
            if(d<0)
            {
                for (it=a[i].begin() ; it != a[i].end() ; it++)
                    s+=fabs(it->second);
                d/=s;
                if(d<min)
                {
                    min=d;
                    l=i;
                }
            }
        }
    }
    return(l);
}

int Simplex::ads_entrant(int n,int l)
{
    int j,k;
    double rap,min;
    min=1e308;
    k=0;
    for(j=1;j<=n;j++)
    {
        if(a[l][j]<0)
        {
            rap=a[0][j]/a[l][j];
            if(rap<min)
            {
                min=rap;
                k=j;
            }
        }
    }
    return(k);
}

void  Simplex::pivotage(int m,int n,int l,int k)
{
    int i,j;
    double pivot,coef;
    pivot=a[l][k];
    for(i=0;i<=m;i++)
        if(i!=l)
        {
            coef=a[i][k]/pivot;
            a[i][k]=-coef;
            for(j=0;j<=n;j++)
                if(j!=k)
                    a[i][j]=a[i][j]-coef*a[l][j];
        }
    coef=1/pivot;
    a[l][k]=coef;
    for(j=0;j<=n;j++)
        if(j!=k)
            a[l][j]=coef*a[l][j];
    i=db[l];
    db[l]=hb[k];
    hb[k]=i;
}

int  Simplex::simplexe_dual(int ineq1, int ineq2, int eq, int n)
{
    try
    {
        int i,j,k,l,phase,m,m1;
        double max;

        m=ineq1+ineq2+eq;

        for(i=ineq1+1;i<=ineq1+ineq2;i++)
            for (it=a[i].begin() ; it != a[i].end() ; it++)
                it->second = -(it->second);
        for(i=1;i<=ineq1+ineq2;i++)
            db[i]=n+i;
        for(i=ineq1+ineq2+1;i<=m;i++)
            db[i]=0;
        for(j=1;j<=n;j++)
            hb[j]=j;
        if(eq!=0)
        {
            for(i=ineq1+ineq2+1;i<=m;i++)
            {
                l=i;
                k=0;
                for(j=1;j<=n;j++)
                    if(a[i][j]!=0)
                        k=j;
                if(k==0)
                {
                    if(a[i][0]!=0)
                        return(2);
                }
                else
                {
                    pivotage(m,n,l,k);
                    hb[k]=hb[n];
                    for(j=0;j<=m;j++)
                        a[j][k]=a[j][n];
                    n-=1;
                }
            }
        }
        m1=m;
        phase=2;
        k=0;
        max=0;
        for(j=1;j<=n;j++)
            if(a[0][j]>max)
            {
                max=a[0][j];
                k=j;
            }
        if(k!=0)
            phase=1;
        l=1;
        if(phase==1)
        {
            m1=m+1;
            for(j=1;j<=n;j++)
                if(a[0][j]>0)
                    a[m1][j]=1;
                else
                    a[m1][j]=0;
            a[m1][0]=0;
            db[m1]=n+m+1;
            pivotage(m1,n,m1,k);
        }
        while(phase<=2)
        {
            do
            {
                l=ads_sortant(m1,n,phase);
                if(l!=0)
                {
                    k=ads_entrant(n,l);
                    if(k==0)
                    {
                        return(1);
                    }
                    pivotage(m1,n,l,k);
                }
            }
            while(l!=0);
            if(phase==1)
            {
                k=0;
                for(j=1;j<=n;j++)
                    if(hb[j]==n+1+m)
                        k=j;
                if(k!=0)
                {
                    if(fabs(a[0][k])>1e-15)
                        return(2);
                    else
                    {
                        for(i=1;i<=m1;i++)
                            if(a[i][k]!=0)
                                l=i;
                        pivotage(m1,n,l,k);
                    }
                }
            }
            if(phase==1)
            {
                for(i=1;i<=m1;i++)
                    if(db[i]==n+m1)
                        l=i;
                db[l]=db[m1];
                for(j=0;j<=n;j++)
                    a[l][j]=a[m1][j];
            }
            phase+=1;
            m1-=1;
        }
        for(i=1;i<=m+n;i++)
            sol[i]=0;
        for(i=1;i<=m;i++)
            sol[db[i]]=a[i][0];
        return (0);
    }
    catch (int erreur)
    {
        if (PFA_VERBOSE)
        {
            cerr << "ERREUR : Simplex::simplexe_dual" << endl;
            switch (erreur)
            {
            default :
                cerr << "Erreur n°" << erreur << endl;
            }
        }
        return erreur;
    }
}

int Simplex::affiche(int ineq1, int ineq2, int nbeq, int nbvar) {
	int i, j;
	
	cout << "--------------------------------" << endl;
	cout << "Affichage du systme :" << endl;
	cout << "--------------------------------" << endl;
	char lettre = 'X';
	cout << " minimize : ";
	for (i=1 ; i<=nbvar ; i++) {
		cout << a[0][i] << lettre++ << " ";
	}
	cout << "::" << endl;
	for (i=1 ; i<=ineq1 ; i++){
		lettre = 'X';
		for (j=1 ; j<=nbvar ; j++) {
			if (a[i][j]>0)
				cout << "+";
			cout << a[i][j] << lettre++ << " ";
		}
		cout << "<=" << a[i][0] << endl;
	}
	for (i=ineq1 + 1; i<= ineq1 + ineq2 ; i++) {
		lettre = 'X';
		for (j=1 ; j<=nbvar ; j++) {
			if (a[i][j]>0)
				cout << "+";
			cout << a[i][j] << lettre++ << " ";
		}
		cout << ">=" << a[i][0] << endl;
	}
	for (i=ineq1 + ineq2 + 1 ; i<=ineq1 + ineq2 + nbeq ; i++) {
		lettre = 'X';
		for (j=1 ; j<=nbvar ; j++) {
			if (a[i][j]>0)
				cout << "+";
			cout << a[i][j] << lettre++ << " ";
		}
		cout << "= " << a[i][0] << endl;
	}
	cout << "--------------------------------" << endl;
	
}

int  Simplex::solve(int ineq1, int ineq2, int eq, int n, bool debug) {
    try {
		const int ntmp = n;
		const int BORNEMAX = 1024;
		int i,j,k;
		int m=ineq1+ineq2+eq;
		int res ;
		
		// MODIFICATION POUR GERER DES VARIABLES NEGATIVES
		for (i=0 ; i <= m ; i++) {
			for (j=1 ; j<=n ; j++) {	
				a[i][j+n] = -a[i][j];
				// cout << "a[" << i <<"][" << j << "]=" << a[i][j];
				// cout << "a[" << i << "][" << j+n << "]=" << a[i][j+n];
			}
			// cout << endl;
		}
		n *= 2;
		// FIN DE LA MODIFICATION
		
		// MODIFICATION POUR EVITER UN BUG DE PRECISION ; on borne
		// toutes les variables
		k = ineq1 + ineq2 + 1;
		for (i=0 ; i<eq ; i++) {
			for (j=0 ; j<=n ; j++) {
				a[k+i+n]=a[k+i];
			}
		}
		for (i=0 ; i<n ; i++) {
			for (j=1 ; j<=n ; j++) {
				a[k+i][j]=0;
			}
			a[k+i][0]=-BORNEMAX;
			a[k+i][i+1]=-1;
		}
		ineq2 = ineq2 + n;
		// FIN DE LA MODIFICATION ANTIBUG DE PRECISION.
	
		if (debug) {
			affiche(ineq1,ineq2,eq,n);
		}

		res = simplexe_dual(ineq1,ineq2,eq,n);
		if (res != 1) {		
			for (i=1 ; i<=ntmp ; i++)
				sol[i] -= sol[i+ntmp];
		
			// // --------- AFfICHAGE DEBUG SOLUTION -------
			// char lettre = 'X';
			// cout << "\nn=" << n << ", m=" << m << endl;
			// cout << "\nSOLUTION : ";
			// for (i=1 ; i<=n ; i++) {
			// 	cout << lettre++ << "=" << sol[i]<< ", ";
			// }
			// cout << endl;
			// // --------- FIN AFFICHAGE DEBUG SOLUTION ------
		}
		return (res);
	}
	catch (int erreur) {
		if (PFA_VERBOSE)
		{
			cerr << "ERREUR : Simplex::solve" << endl;
			switch (erreur)
			{
				default :
					cerr << "Erreur n°" << erreur << endl;
			}
		}
		return erreur;
	}
}


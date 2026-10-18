/* 
   Prototypes de fonctions
*/
int TraitementCoordonnees(Tcl_Interp *interp,
			  Tcl_Obj **OLinePtr,
			  int nbIndividus,
			  int nbDimensions,
			  individu_t ***individus, 
			  int typeDonnees);

int TraitementDistances(Tcl_Interp *interp,
			Tcl_Obj **OLinePtr,
			int nbIndividus,
			int nbDimensions,
			individu_t ***individus,
			double ***distances);

int TraitementAlignement(Tcl_Interp *interp,
			 Tcl_Obj **OLine,
			 int nbIndividus,
			 int *longueurAlignement,
			 individu_t ***individus,
			 char ***sequences);

/* fin prototypage */

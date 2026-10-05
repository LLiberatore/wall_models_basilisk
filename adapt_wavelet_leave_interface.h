/**
# adapt_wavelet_leave_interface_limited()

Variante di `adapt_wavelet_leave_interface()` in cui il livello massimo di
raffinamento NON è un singolo intero globale, ma dipende dalla posizione
nel dominio, tramite una funzione definita dall'utente:

  int local_maxlevel (double x, double y, double z);

che va dichiarata/definita PRIMA di includere questo header (in 2D il
parametro z esiste comunque, vale 0 per convenzione Basilisk, e può
essere semplicemente ignorato nel corpo della funzione).

Rispetto a `adapt_wavelet_leave_interface`:

- il criterio wavelet confronta `level` con `local_maxlevel(x,y,z)`
  invece che con un `p.maxlevel` fisso;
- le celle di interfaccia (0 < vol_frac < 1, con eventuale padding di
  vicini) continuano a essere SOLO forzate verso il refine finche' non
  raggiungono il loro maxlevel locale — esattamente come nell'originale
  `adapt_wavelet_leave_interface`. Una volta raggiunto il maxlevel della
  propria zona, la cella di interfaccia non viene piu' toccata (ne'
  refine ne' coarsening) da questo blocco. Non c'e' bisogno di gestire
  un "cambio di zona" a runtime, perche' la zona di una cella dipende
  solo dalla sua posizione (x,y), che e' fissa: una cella non puo' mai
  "entrare" in un'altra zona durante la simulazione.
- il bilanciamento 2:1 tra zone a maxlevel diverso è garantito
  nativamente dalla struttura ad albero di Basilisk (refine_cell /
  coarsen_cell + mpi_boundary_refine/coarsen), non richiede logica
  aggiuntiva finché il salto di maxlevel tra zone confinanti è <= 1.
*/

struct Adapt_leave_interface_limited {
  scalar * slist;     // scalari usati per la stima wavelet (es. {u,f})
  scalar * vol_frac;  // scalari di tipo volume-fraction/cs che segnano
                       // le celle di interfaccia da preservare (es. {cs})
  double * max;        // tolleranza per ciascuno scalare di slist
  int minlevel;         // livello minimo di raffinamento (default 1)
  int padding;           // n. di celle vicine da includere attorno
                          // all'interfaccia (default 0)
  scalar * list;          // campi da aggiornare nell'albero (default all)
};

astats adapt_wavelet_leave_interface_limited (struct Adapt_leave_interface_limited p)
{
  scalar * list = p.list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (p.slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
        list = list_add (list, s);
    }
    boundary (list);
    scalar * listr = list_concat (p.slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    if (!is_constant(s) && s.restriction != no_restriction)
      listc = list_add (listc, s);

  if (p.minlevel < 1)
    p.minlevel = 1;

  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);

  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);

      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++;
        }
        continue;
      }
      else { // !is_leaf (cell)
        if (cell.flags & refined) {
          // cella già raffinata in questo passo, salta i figli
          cell.flags &= ~too_coarse;
          continue;
        }

        // controlla se la cella o uno dei figli è locale (MPI)
        bool local = is_local(cell);
        if (!local)
          foreach_child()
            if (is_local(cell)) {
              local = true; break;
            }

        if (local) {

          // --- 1. criterio wavelet standard, con maxlevel per-cella ---
          int i = 0;
          static const int just_fine = 1 << (user + 3);
          for (scalar s in p.slist) {
            double max = p.max[i++], sc[1 << dimension];
            int c = 0;
            foreach_child()
              sc[c++] = s[];
            s.prolongation (point, s);
            c = 0;
            foreach_child() {
              int cellMAX = local_maxlevel (x, y, z);
              double e = fabs(sc[c] - s[]);
              if (e > max && level < cellMAX) {
                cell.flags &= ~too_fine;
                cell.flags |= too_coarse;
              }
              else if ((e <= max/1.5 || level > cellMAX) &&
                       !(cell.flags & (too_coarse|just_fine))) {
                if (level >= p.minlevel)
                  cell.flags |= too_fine;
              }
              else if (!(cell.flags & too_coarse)) {
                cell.flags &= ~too_fine;
                cell.flags |= just_fine;
              }
              s[] = sc[c++];
            }
          }

          // --- 2. vincolo interfaccia (vol_frac), con maxlevel per-cella.
          //        Come nell'originale: si forza SOLO il refine verso il
          //        maxlevel locale finche' non lo si raggiunge; una volta
          //        raggiunto, la cella non viene piu' toccata (ne' refine
          //        ne' coarsening) da questo blocco. ---
          foreach_child() {
            for (scalar vf in p.vol_frac) {
              int cellMAX = local_maxlevel (x, y, z);
              if (vf[] > 1e-4 && vf[] < 1. - 1e-4 && level < cellMAX) {
                cell.flags |= too_coarse;
                cell.flags &= ~too_fine;
                cell.flags &= ~just_fine;
                if (p.padding > 0)
                  foreach_neighbor (p.padding) {
                    // il maxlevel va ricalcolato nelle coordinate del
                    // vicino: se il padding sconfina in un'altra zona,
                    // il vicino viene forzato al SUO maxlevel locale,
                    // non a quello della cella di interfaccia centrale
                    int nMAX = local_maxlevel (x, y, z);
                    if (level < nMAX) {
                      cell.flags |= too_coarse;
                      cell.flags &= ~too_fine;
                      cell.flags &= ~just_fine;
                    }
                  }
              }
            }
          }

          foreach_child() {
            int cellMAX = local_maxlevel (x, y, z);
            cell.flags &= ~just_fine;
            if (!is_leaf(cell)) {
              cell.flags &= ~too_coarse;
              if (level >= cellMAX)
                cell.flags |= too_fine;
            }
            else if (!is_active(cell))
              cell.flags &= ~too_coarse;
          }
        }
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);

  for (int l = depth(); l >= 0; l--) {   // coarsening 
    foreach_cell()
      if (!is_boundary(cell)) {
        if (level == l) {
          if (!is_leaf(cell)) {
            if (cell.flags & refined)
              cell.flags &= ~(refined|too_fine);
            else if (cell.flags & too_fine) {
              if (is_local(cell) && coarsen_cell (point, listc))
                st.nc++;
              cell.flags &= ~too_fine; // non coarsenare il genitore
            }
          }
          if (cell.flags & too_fine)
            cell.flags &= ~too_fine;
          else if (level > 0 && (aparent(0).flags & too_fine))
            aparent(0).flags &= ~too_fine;
          continue;
        }
        else if (is_leaf(cell))
          continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != p.list)
    free (list);

  return st;
}

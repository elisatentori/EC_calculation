/*=============================================================================
Code update
Copyright (c) 2024 – University of Padua, Padova Neuroscience Center

Author: Elisa Tentori (elisa.tentori@phd.unipd.it)

Code modified in order to calculate TE and to count delayed coincidences between 
pre-synaptic jittered timeseries and post-synaptic not-jittered timeseries.
This C code is called by ASDFTE_perm function (should be in the same folder).

To compile in matlab: mex transent_CN_perm.c


NOTES:
(1) In this code x=i is the POST synaptic neuron series of spike-times,
    While y=j is the PRE synaptic neuron series of spike-times
(2) The original transent_1(..) function calculates the TE as
    TE(x(t),x(t-1),y(t-d)). The choice, motivated by the autors in Ito et al.-2011,
    is based on the assumption that the only past of the post-synaptic neuron
    that counts is the one related to its last 1 ms time window.
    Here the function has been modified to calculate TE as
    TE(x(t),x(t-d),y(t-d)). Indeed, we prefer to evaluate also x at time d.
=============================================================================*/

/*=============================================================================
Original code:
Copyright (c) 2011, The Trustees of Indiana University
All rights reserved.

Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu)
=============================================================================*/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <memory.h>


typedef int TimeType;

/* TE order 1 */
void transent_1
(const mxArray *all_series_y, // y : pre-synaptic neuron
 const mxArray *all_series_x, // x : post-synaptic neuron
 const mwSize series_count,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result, double *coincidence_result);   //# Modifica: Aggiunto coincidence_result come parametro //#


/* TE higher order than 1 */
void transent_ho
(const mxArray *all_series_y,
 const mxArray *all_series_x,
 const mwSize series_count,
 const unsigned int x_order, const unsigned int y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result, double *coincidence_result);   // Modifica: Aggiunto coincidence_result come parametro //#




/* Interface betw C++ and matlab */
/* transent_perm(asdf->jittered, asdf2, y_delay, x_order=1, y_order=1) */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

  mwSize series_count;
  mwSize series_count2;
  unsigned int x_order, y_order;
  TimeType y_delay, duration, duration2;
  double *data_ptr;
  double *data_ptr2;
  mxArray *array_ptr;
  mxArray *array_ptr2;

  if (nlhs > 2) {  // Modifica: Cambiato il controllo per due output //#
		mexErrMsgTxt("Expected only two output arguments");  //#
  }

  /* Arguments: all_series_y, all_series_x, y_delay */
  if (nrhs < 3) {
		mexErrMsgTxt("Expected at least three input arguments");
  }

  if (!mxIsCell(prhs[0])) {
		mexErrMsgTxt("First argument must be a cell array");
  }
  
  if (!mxIsCell(prhs[1])) {
		mexErrMsgTxt("Second argument must be a cell array");
  }

  /* --- Extract arguments --- */
  
  /* Number of timeseries per Cell array */
  series_count = mxGetNumberOfElements(prhs[0]) - 2; /* Last two cells hold info */
  // series_count gives the number of elements (rows) in the specified mxArray.
  series_count2 = mxGetNumberOfElements(prhs[1]) - 2;
  
  if (series_count!=series_count2) {
		mexErrMsgTxt("Second argument must have same number of series than the first argument");
  }


  /* Extract duration of first time-series arrays */
  array_ptr = mxGetCell(prhs[0], series_count + 1);
  // mwPointer mxGetCell(pm, index)
  // pm: Pointer to a cell mxArray  –  index: Number of elements in the cell mxArray between
  //                                          the first element and the desired one.
  data_ptr = mxGetPr(array_ptr);
  // input: Pointer to a MATLAB array of type mxDOUBLE_CLASS, specified as mxArray *.
  // output: Pointer to the data array within an mxArray, specified as mxDouble *.
  duration = (TimeType)data_ptr[1];

  /* Extract duration of second time-series arrays */
  array_ptr2 = mxGetCell(prhs[1], series_count2 + 1);
  data_ptr2 = mxGetPr(array_ptr2);
  duration2 = (TimeType)data_ptr2[0];

  /* Check: duration has to be the same in the 2 timeseries
  if (duration != duration2) {
		mexErrMsgTxt("Second argument must have same dimension than the first argument");
  }*/
  
  
  /* Extract y_time_delay */
  data_ptr = mxGetPr(prhs[2]); /* points to the third argument, that is y_delay */
  y_delay = (TimeType)data_ptr[0];
  
  
  /* x_order and y_order */
  /*
     Nota che se inserisci solo il quarto argomento e non anche il quinto, risetta sia il quarto che
     il quinto (cioè x_order e y_order) a 1 di default. Quindi se li specifichi, li devi specificare
     entrambi, altrimenti se ne inizializzi solo uno, te lo rimette a 1 di default insieme all'altro.
  */
  if (nrhs == 5) {
    data_ptr = mxGetPr(prhs[3]);
    x_order = (unsigned int)data_ptr[0];

    data_ptr = mxGetPr(prhs[4]);
    y_order = (unsigned int)data_ptr[0];
  }
  else {
    x_order = 1; // default values
    y_order = 1;
  }

  /* Create result matrix */
  plhs[0] = mxCreateDoubleMatrix(series_count, series_count, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(series_count, series_count, mxREAL);  //# Modifica: Creazione matrice coincidence_result

  /* Do calculation */
  data_ptr = mxGetPr(plhs[0]); /* Output object */
  double *coincidence_result = mxGetPr(plhs[1]);  //# Modifica: Puntatore a coincidence_result

  if ((x_order == 1) && (y_order == 1)) {
    transent_1(prhs[0], prhs[1], series_count,
               y_delay, duration,
               data_ptr, coincidence_result);  // Modifica: Passaggio di coincidence_result //#
  } else {
    transent_ho(prhs[0], prhs[1], series_count,
                x_order, y_order,
                y_delay, duration,
                data_ptr, coincidence_result);  // Modifica: Passaggio di coincidence_result //#
  }
}


/* *************************************** */

double Log2(double n) {
	return log(n) / log(2.0);
}

/* *************************************** */

/* Computes the first-order transfer entropy matrix for all pairs. */
void transent_1
(const mxArray *all_series_y, const mxArray *all_series_x,
 const mwSize series_count,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result, double *coincidence_result) {  // Modifica: Aggiunto coincidence_result come parametro //#


  /* Constants */
  const unsigned int x_order = 1, y_order = 1,                
               num_series = 3,
               num_counts = 8,
               num_x = 4,
               num_y = 2;

  /* Locals */
  TimeType counts[8];
  unsigned long code;
  long k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  double *ord_iter[3];
  double *ord_end[3];

  TimeType ord_times[3];
  TimeType ord_shift[3];

  const unsigned int window = y_order + y_delay;
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  mxArray *array_ptr;
  double *i_series, *j_series;
  mwSize i_size, j_size;
  mwIndex i, j;

  /* MATLAB is column major */
  for (j = 0; j < series_count; ++j) {
    for (i = 0; i < series_count; ++i) {

      /* Extract series */
      array_ptr = mxGetCell(all_series_x, i);
      i_size = mxGetNumberOfElements(array_ptr);
      i_series = mxGetPr(array_ptr);
    
      /* all_series_y is the jittered pre-synaptic timeserie */
      array_ptr = mxGetCell(all_series_y, j);
      j_size = mxGetNumberOfElements(array_ptr);
      j_series = mxGetPr(array_ptr);

      if ((i_size == 0) || (j_size == 0)) {
        te_result[(i * series_count) + j] = 0;
        coincidence_result[(i * series_count) + j] = 0;  //# Modifica: Inizializza la cella delle coincidenze a zero
        continue;
      }

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) { // ORIGINAL [default]
      //for (k = 0; k < (x_order); ++k) {  // MODIFIED
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = -k;
        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);
      TimeType coincidence_count = 0;  //# Modifica: Inizializza il conteggio delle coincidenze


      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {      
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = (TimeType)*(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        
        //# Coincidence between X(t) and Y(t-y_delay) 
        if (code == 0b101 || code == 0b111) {//if (code == 0b101) { // alternative
            ++coincidence_count;
        }

        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * Log2(prob_2 / prob_3));
      }

      /* MATLAB is column major, but flipped for compatibility */
      te_result[(i * series_count) + j] = te_final;
      coincidence_result[(i * series_count) + j] = coincidence_count;  //# Modifica: Assegna il conteggio delle coincidenze


    } /* for i */

  } /* for j */
 
} /* transent_1 */


/* Computes the higher-order transfer entropy matrix for all pairs. */
void transent_ho
(const mxArray *all_series_y, const mxArray *all_series_x,
 const mwSize series_count,
 const unsigned int x_order, const unsigned int y_order,
 const TimeType y_delay,
 const TimeType duration,
 double *te_result, double *coincidence_result) {  //# Modifica: Aggiunto coincidence_result come parametro


  /* Constants */
  const unsigned int num_series = 1 + y_order + x_order,
               num_counts = (unsigned int)pow(2, num_series),
               num_x = (unsigned int)pow(2, x_order + 1),
               num_y = (unsigned int)pow(2, y_order);

  /* Locals */
  TimeType *counts = (TimeType*)malloc(sizeof(TimeType) * num_counts);
  unsigned long code;
  long k, l, idx, c1, c2;
  double te_final, prob_1, prob_2, prob_3;

  double **ord_iter = (double**)malloc(sizeof(double*) * num_series);
  double **ord_end = (double**)malloc(sizeof(double*) * num_series);

  TimeType *ord_times = malloc(sizeof(TimeType) * num_series);
  TimeType *ord_shift = malloc(sizeof(TimeType) * num_series);

  const unsigned int window = (y_order + y_delay) > (x_order + 1) ? (y_order + y_delay) : (x_order + 1);
  const TimeType end_time = duration - window + 1;
  TimeType cur_time, next_time;

  /* Calculate TE */
  mxArray *array_ptr;
  double *i_series, *j_series; /* i = pre-synaptic; j = post-synaptic*/
  mwSize i_size, j_size;
  mwIndex i, j;

  /* MATLAB is column major */
  for (j = 0; j < series_count; ++j) {
    for (i = 0; i < series_count; ++i) {

      /* Extract series */
      array_ptr = mxGetCell(all_series_x, i);
      i_size = mxGetNumberOfElements(array_ptr);
      i_series = mxGetPr(array_ptr);

      array_ptr = mxGetCell(all_series_y, j);
      j_size = mxGetNumberOfElements(array_ptr);
      j_series = mxGetPr(array_ptr);

      if ((i_size == 0) || (j_size == 0)) {
        te_result[(i * series_count) + j] = 0;
        coincidence_result[(i * series_count) + j] = 0;  //# Modifica: Inizializza la cella delle coincidenze a zero
		continue;
      }

      /* Order is x^(k+1), y^(l) */
      idx = 0;

      /* x^(k+1) */
      for (k = 0; k < (x_order + 1); ++k) {
        ord_iter[idx] = i_series;
        ord_end[idx] = i_series + i_size;
        ord_shift[idx] = (window - 1) - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* y^(l) */
      for (k = 0; k < y_order; ++k) {
        ord_iter[idx] = j_series;
        ord_end[idx] = j_series + j_size;
        ord_shift[idx] = (window - 1) - y_delay - k;

        while ((TimeType)*(ord_iter[idx]) < ord_shift[idx] + 1) {
          ++(ord_iter[idx]);
        }

        ord_times[idx] = (TimeType)*(ord_iter[idx]) - ord_shift[idx];
        ++idx;
      }

      /* Count spikes */
      memset(counts, 0, sizeof(TimeType) * num_counts);
      TimeType coincidence_count = 0;  //# Modifica: Inizializza il conteggio delle coincidenze

      /* Get minimum next time bin */
      cur_time = ord_times[0];
      for (k = 1; k < num_series; ++k) {
        if (ord_times[k] < cur_time) {
          cur_time = ord_times[k];
        }
      }

      while (cur_time <= end_time) {

        code = 0;
        next_time = end_time + 1;

        /* Calculate hash code for this time bin */
        for (k = 0; k < num_series; ++k) {
          if (ord_times[k] == cur_time) {       
            code |= 1 << k;

            /* Next spike for this neuron */
            ++(ord_iter[k]);

            if (ord_iter[k] == ord_end[k]) {
              ord_times[k] = end_time + 1;
            }
            else {
              ord_times[k] = (TimeType)*(ord_iter[k]) - ord_shift[k];
            }
          }

          /* Find minimum next time bin */
          if (ord_times[k] < next_time) {
            next_time = ord_times[k];
          }
        }

        ++(counts[code]);
        
        //# Aggiunto per conta coincidenze
        int y_delay_bit = x_order + y_delay;
        if ((code & (1 << 0)) && (code & (1 << y_delay_bit))) {
            ++coincidence_count;
        }
          
        cur_time = next_time;

      } /* while spikes left */

      /* Fill in zero count */
      counts[0] = end_time;
      for (k = 1; k < num_counts; ++k) {
        counts[0] -= counts[k];
      }

      /* ===================================================================== */

      /* Use counts to calculate TE */
      te_final = 0;

      /* Order is x^(k), y^(l), x(n+1) */
      for (k = 0; k < num_counts; ++k) {
        prob_1 = (double)counts[k] / (double)end_time;

        if (prob_1 == 0) {
          continue;
        }

        prob_2 = (double)counts[k] / (double)(counts[k] + counts[k ^ 1]);

        c1 = 0;
        c2 = 0;

        for (l = 0; l < num_y; ++l) {
          idx = (k & (num_x - 1)) + (l << (x_order + 1));
          c1 += counts[idx];
          c2 += (counts[idx] + counts[idx ^ 1]);
        }

        prob_3 = (double)c1 / (double)c2;

        te_final += (prob_1 * Log2(prob_2 / prob_3));
      }

      /* MATLAB is column major, but flipped for compatibility */
      te_result[(i * series_count) + j] = te_final;
      coincidence_result[(i * series_count) + j] = coincidence_count;  // Modifica: Assegna il conteggio delle coincidenze //#


    } /* for i */

  } /* for j */

  /* Clean up */
  free(counts);
  free(ord_iter);
  free(ord_end);
  free(ord_times);
  free(ord_shift);
 
} /* transent_ho */


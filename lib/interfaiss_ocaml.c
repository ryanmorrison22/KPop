/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include "interfaiss.h"
#define CAML_NAME_SPACE
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>
#include <caml/fail.h>
#include <caml/bigarray.h>

CAMLprim value InterfaissCreate(value o_index_type,value o_n) {
  CAMLparam2(o_index_type,o_n);
  CAMLlocal1(o_res);
  /* Remember that the OCaml definition of index_t is
    type index_t =
      | Flat
      | PQ of int * int
      | HNSW of int * int */
  index_t* res;
  if (o_index_type==Val_int(0))
    res=interfaiss_create_flat_index(Int_val(o_n));
  else if (Tag_val(o_index_type)==0)
    res=interfaiss_create_PQ_index(Int_val(o_n),Int_val(Field(o_index_type,0)),Int_val(Field(o_index_type,1)));
  else if (Tag_val(o_index_type)==1)
    res=interfaiss_create_HNSW_index(Int_val(o_n),Int_val(Field(o_index_type,0)),Int_val(Field(o_index_type,1)));
  else
    caml_failwith("InterfaissMake");
  o_res=caml_alloc(1,Abstract_tag);
  *((index_t**)Data_abstract_val(o_res))=res;
  CAMLreturn(o_res); 
}

CAMLprim value InterfaissAdd(value o_idx,value o_vectors) {
  CAMLparam2(o_idx,o_vectors);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  interfaiss_add_data_to_index(idx,d,n,Caml_ba_data_val(o_vectors));
  return Val_int(0);
}

CAMLprim value InterfaissTrain(value o_idx,value o_vectors) {
  CAMLparam2(o_idx,o_vectors);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  interfaiss_train_index(idx,d,n,Caml_ba_data_val(o_vectors));
  return Val_int(0);
}

CAMLprim value InterfaissQuery(value o_idx,value o_vectors,value o_k) {
  CAMLparam3(o_idx,o_vectors,o_k);
  CAMLlocal1(o_res);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  int k=Int_val(o_k);
  if (k<0)
    caml_failwith("InterfaissQuery");

  interfaiss_query_index(idx,d,n,Caml_ba_data_val(o_vectors),k, float* distances, idx_t* indices)


  CAMLreturn(o_res);
}

CAMLprim value InterfaissDelete(value o_idx) {
  CAMLparam1(o_idx);
  interfaiss_free_index(*((index_t**)Data_abstract_val(o_idx)));
  return Val_int(0);
}

/*
CAMLprim value Bundl_Aligner(
    value o_sequences,value aligner_type,value num_iterations,
    value flag_robust_alignment,value flag_time,value flag_verbose,value num_threads) {
  CAMLparam5(o_sequences,aligner_type,num_iterations,flag_robust_alignment,flag_time);
  CAMLxparam2(flag_verbose,num_threads);
  CAMLlocal1(o_aligned);
  const size_t n=(size_t)Wosize_val(o_sequences);
  char const** c_sequences=calloc(n,sizeof(char*));
  size_t i;
  for (i=0;i<n;++i)
    c_sequences[i]=String_val(Field(o_sequences,i));
  char** c_aligned=
    bundl_run(c_sequences,(unsigned int)n,(unsigned int)Int_val(num_threads),
              (bundl_aligner_index_t)Int_val(aligner_type),(unsigned int)Int_val(num_iterations),
              Bool_val(flag_time),Bool_val(flag_verbose),Bool_val(flag_robust_alignment));
  if (c_aligned==NULL)
    o_aligned=Val_int(0);
  else {
    o_aligned=caml_copy_string_array((char const**)c_aligned);
    for (i=0;i<n;++i)
      free(c_aligned[i]);
    free(c_aligned);
  }
  free(c_sequences);
  CAMLreturn(o_aligned);
}
CAMLprim value Bundl_Aligner_Bytecode(value* argv,int argn) {
  return Bundl_Aligner(argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
}
*/


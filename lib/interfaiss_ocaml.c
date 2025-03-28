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
      | HNSW of int */
  index_t* res;
  if (o_index_type==Val_int(0))
    res=interfaiss_create_flat_index(Int_val(o_n));
  else if (Tag_val(o_index_type)==0)
    res=interfaiss_create_PQ_index(Int_val(o_n),Int_val(Field(o_index_type,0)),Int_val(Field(o_index_type,1)));
  else if (Tag_val(o_index_type)==1)
    res=interfaiss_create_HNSW_index(Int_val(o_n),Int_val(Field(o_index_type,0)),0);
  else
    caml_failwith("InterfaissCreate");
  o_res=caml_alloc(1,Abstract_tag);
  *((index_t**)Data_abstract_val(o_res))=res;
  CAMLreturn(o_res);
}

CAMLprim value InterfaissAdd(value o_idx,value o_vectors) {
  CAMLparam2(o_idx,o_vectors);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  interfaiss_add_data_to_index(idx,d,(idx_t)n,Caml_ba_data_val(o_vectors));
  CAMLreturn(Val_unit);
}

CAMLprim value InterfaissTrain(value o_idx,value o_vectors) {
  CAMLparam2(o_idx,o_vectors);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  interfaiss_train_index(idx,d,(idx_t)n,Caml_ba_data_val(o_vectors));
  CAMLreturn(Val_unit);
}

CAMLprim value InterfaissQuery(value o_idx,value o_vectors,value o_k) {
  CAMLparam3(o_idx,o_vectors,o_k);
  CAMLlocal1(o_res);
  index_t* idx=*((index_t**)Data_abstract_val(o_idx));
  int n=Caml_ba_array_val(o_vectors)->dim[0];
  int d=Caml_ba_array_val(o_vectors)->dim[1];
  idx_t k=(idx_t)Int_val(o_k);
  if (k<=0)
    caml_failwith("InterfaissQuery");
  float* distances=NULL;
  idx_t* indices=NULL;
  interfaiss_query_index(idx,d,(idx_t)n,Caml_ba_data_val(o_vectors),&k,&distances,&indices);
  o_res=caml_alloc_tuple(2);
  long dims[2];
  dims[0]=(long)n;
  dims[1]=(long)k;
  // BEWARE: THIS STATEMENT DEPENDS ON THE TYPE USED BY faiss
  Store_field(o_res,0,caml_ba_alloc(CAML_BA_INT64|CAML_BA_C_LAYOUT,2,indices,dims));
  Store_field(o_res,1,caml_ba_alloc(CAML_BA_FLOAT32|CAML_BA_C_LAYOUT,2,distances,dims));
  CAMLreturn(o_res);
}

CAMLprim value InterfaissDelete(value o_idx) {
  CAMLparam1(o_idx);
  interfaiss_free_index(*((index_t**)Data_abstract_val(o_idx)));
  CAMLreturn(Val_unit);
}


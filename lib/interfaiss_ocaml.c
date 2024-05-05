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
//#include <caml/fail.h>

CAMLprim value InterfaissMake(value o_index_type,value o_vector) {





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


/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/experimental/emap_lite.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/utils/struct_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/write_dot.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/verilog_reader.hpp>

#include <experiments.hpp>

std::string const mcnc_library = "GATE   inv1    1  O=!a;             PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2  O=!a;             PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   inv3    3  O=!a;             PIN * INV 3 999 1.1 0.09 1.1 0.09\n"
                                 "GATE   inv4    4  O=!a;             PIN * INV 4 999 1.2 0.07 1.2 0.07\n"
                                 "GATE   nand2   2  O=!(a*b);         PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   nand3   3  O=!(a*b*c);       PIN * INV 1 999 1.1 0.3 1.1 0.3\n"
                                 "GATE   nand4   4  O=!(a*b*c*d);     PIN * INV 1 999 1.4 0.4 1.4 0.4\n"
                                 "GATE   nor2    2  O=!(a+b);         PIN * INV 1 999 1.4 0.5 1.4 0.5\n"
                                 "GATE   nor3    3  O=!(a+b+c);       PIN * INV 1 999 2.4 0.7 2.4 0.7\n"
                                 "GATE   nor4    4  O=!(a+b+c+d);     PIN * INV 1 999 3.8 1.0 3.8 1.0\n"
                                 "GATE   and2    3  O=a*b;            PIN * NONINV 1 999 1.9 0.3 1.9 0.3\n"
                                 "GATE   or2     3  O=a+b;            PIN * NONINV 1 999 2.4 0.3 2.4 0.3\n"
                                 "GATE   xor2a   5  O=a*!b+!a*b;      PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "#GATE  xor2b   5  O=!(a*b+!a*!b);   PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   xnor2a  5  O=a*b+!a*!b;      PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                 "#GATE  xnor2b  5  O=!(a*!b+!a*b);   PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                                 "GATE   aoi21   3  O=!(a*b+c);       PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                 "GATE   aoi22   4  O=!(a*b+c*d);     PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "GATE   oai21   3  O=!((a+b)*c);     PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                                 "GATE   oai22   4  O=!((a+b)*(c+d)); PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "GATE   or6     6  O=(a+b+c+d+e+f);  PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "#GATE   and6    6  O=(a*b*c*d*e*f);  PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                                 "GATE   buf     2  O=a;              PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0  O=CONST0;\n"
                                 "GATE   one     0  O=CONST1;";

using namespace experiments;
using namespace mockturtle;

std::string bench_name = "custom_benchmark";
std::string bench_full_name = "/home/radi/RA/mockturtle/experiments/benchmarks/" + bench_name;
std::string res_full_name = "result.v";
std::string dot_aig_full_name = "aig.dot";
std::string dot_res_full_name = "result.dot";




int main()
{
  
  experiment<std::string, uint32_t, double, uint32_t, double, float, bool> exp(
      "emap", "benchmark", "size", "area_after", "depth", "delay_after", "runtime", "cec" );

  fmt::print( "[i] processing technology library\n" );

  /* library to map to technology */
  std::vector<gate> gates;
  //std::ifstream in( "/home/radi/RA/mockturtle/asap_smp.genlib" ); //library or9
  std::stringstream in ( mcnc_library );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

   //compute struct library
  struct_library<6> str_lib( gates );

  tech_library_params tps;
  tps.verbose = true;
  tech_library<6, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : iscas_benchmarks( ) ) //iscas_benchmarks() //or 2nd epfl benchmarks // experiments::bar
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;

    if ( lorina::read_aiger( benchmark_path(benchmark), aiger_reader( aig ) ) != lorina::return_code::success )
     continue;

    /*if ( lorina::read_verilog( bench_full_name + ".v", verilog_reader( aig ) ) != lorina::return_code::success )
      continue;*/
    
    //verification format
    write_aiger(aig, bench_full_name + ".aig" );

    //aig_balance( aig );

    write_dot(aig, dot_aig_full_name);

    write_verilog(aig, "aig_balanced.v");

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    emap_lite_params ps;
    ps.verbose = true;
    emap_lite_stats st;

    //cut size = 6 and Ninputs = 6
    binding_view<klut_network> res = emap_lite<aig_network, 6, 6>( aig, tech_lib, str_lib, ps, &st );

    res.report_gates_usage();

    write_verilog_with_binding(res, res_full_name);

    write_dot(res, dot_res_full_name);

    //bool const cec = benchmark != "hyp" ? abc_cec( res, bench_name ) : true;

    bool const cec = benchmark != "hyp" ? abc_cec( res, benchmark ) : true;

    exp( benchmark, size_before, st.area, depth_before, st.delay, to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}

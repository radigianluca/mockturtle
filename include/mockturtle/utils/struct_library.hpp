
#pragma once

#include <cassert>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <numeric>
#include <bits/stdc++.h>


#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include <parallel_hashmap/phmap.h>

#include "super_utils.hpp"
#include "tech_library.hpp"
#include "../io/genlib_reader.hpp"
#include "../io/super_reader.hpp"

namespace mockturtle
{

/*
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
                                 "GATE   buf     2  O=a;              PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0  O=CONST0;\n"
                                 "GATE   one     0  O=CONST1;";
*/

/*! \brief Library of gates for Boolean matching
 *
 * This class creates a technology library from a set
 * of input gates. Each NP- or P-configuration of the gates
 * are enumerated and inserted in the library.
 * 
 * The configuration is selected using the template
 * parameter `Configuration`. P-configuration is suggested
 * for big libraries with few symmetric gates. The template
 * parameter `NInputs` selects the maximum number of variables
 * allowed for a gate in the library.
 * 
 * The library can be generated also using supergates definitions.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      std::vector<gate> gates;
      lorina::read_genlib( "file.genlib", genlib_reader( gates ) );
      // standard library
      mockturtle::tech_library lib( gates );

      super_lib supergates_spec;
      lorina::read_super( "file.super", super_reader( supergates_spec ) );
      // library with supergates
      mockturtle::tech_library lib_super( gates, supergates_spec );
   \endverbatim
 */
template<unsigned NInputs = 4u>
class struct_library
{

  enum class node_type
  {
  none,
  pi_,
  and_,
  or_,
  mux_,
  xor_
  };

  struct signal
  {
    uint64_t inv : 1;
    uint64_t index :63;
  };

  struct dsd_node //dsd dec
  {
    node_type type;

    int index;

    /*std::vector<bool> inv = {false, false};

    std::vector<dsd_node> fanin = {};*/

    std::vector<signal> fanin = {};
  };

  using supergates_list_t = std::vector<supergate<NInputs>>;
  using tt_hash = kitty::hash<kitty::static_truth_table<NInputs>>;
  using lib_t = phmap::flat_hash_map<kitty::static_truth_table<NInputs>, supergates_list_t, tt_hash>;
  using lib_rule = phmap::flat_hash_map<kitty::dynamic_truth_table, std::vector<dsd_node>, kitty::hash<kitty::dynamic_truth_table>>;

public:

  std::string to_string(node_type t)
  {
    if(t == node_type::and_)
      return "*";
    if(t == node_type::or_)
      return "+";
    if(t == node_type::mux_)
      return "+";
    if(t == node_type::xor_)
      return "xor";
    if(t == node_type::pi_)
      return "pi";
    if(t == node_type::none)
      return "none";
  }

  explicit struct_library( std::vector<gate> const& gates, tech_library_params const ps = {}, super_lib const& supergates_spec = {} )
      : _gates( gates ),
        _supergates_spec( supergates_spec ),
        _ps( ps ),
        _super( _gates, _supergates_spec ),
        _use_supergates( false ),
        _super_lib(),
        _rules_map()
  {
    generate_library();
  }

  explicit struct_library( std::vector<gate> const& gates, super_lib const& supergates_spec, tech_library_params const ps = {} )
      : _gates( gates ),
        _supergates_spec( supergates_spec ),
        _ps( ps ),
        _super( _gates, _supergates_spec, super_utils_params{ps.verbose} ),
        _use_supergates( true ),
        _super_lib(),
        _rules_map()
  {
    generate_library();
  }

  /*! \brief Get the gates matching the function.
   *
   * Returns a list of gates that match the function represented
   * by the truth table.
   */
  const std::vector<kitty::top_decomposition>* get_supergates( kitty::static_truth_table<NInputs> const& tt ) const
  {
    auto match = _super_lib.find( tt );
    if ( match != _super_lib.end() )
      return &match->second;
    return nullptr;
  }

  /*! \brief Get the rules matching the function.
   *
   * Returns a list of gates that match the function represented
   * by the truth table.
   */
  std::vector<dsd_node> get_rules( kitty::dynamic_truth_table const& tt )
  {
    auto match = _rules_map.find( tt );
    if ( match != _rules_map.end() )
      return match->second;
    return {};
  }

  /*! \brief Returns the original gates. */
  const std::vector<gate> get_gates() const
  {
    return _gates;
  }

private:
  void generate_library()
  {

    std::vector<std::vector<dsd_node>> rules = {};
    //std::vector<gate_struct> bttm_rules = {};

    auto supergates = _super.get_super_library();
    uint32_t const standard_gate_size = _super.get_standard_library_size();

    //sort increasing order of area
    sort(supergates.begin(), supergates.end(), 
      []( auto const& gate1, auto const& gate2) -> bool {
      return gate1.area < gate2.area;
    });

    //decomposition
    for(auto& gate : supergates)
    {
      if(gate.num_vars > 1)
      {
        std::vector<dsd_node> rule = {};
        std::vector<int> support = {};
        for(int i = 0; i < gate.num_vars;i++)
        {
          rule.push_back({node_type::pi_, i, {}});
          support.push_back(i);
        }
        auto cpy = gate.function;
        compute_rule(cpy, support, rule);
        _rules_map.insert({gate.function, rule});
        rec_print_dsd(rule, rule[rule.size()-1]);
        std::cout<<"\n";
        rules.push_back(rule);
        //_rules_map.insert( {gate.function, rule} );
      }
    }
    std::cout << "\n";
    
  }

  void print_dsd_node(dsd_node n)
  {
    std::cout << "Type " << to_string(n.type) << " index " << n.index << "\n";
    for(auto elem : n.fanin)
    {
      if(elem.inv)
        std::cout << "!";
      std::cout << elem.index << "\t";
    }
    std::cout << "\n";
  }

  void rec_print_dsd(std::vector<dsd_node> rule, dsd_node n)
  {
    if(n.type == node_type::pi_)
    {
      std::cout<< char('a' + n.index);
      return;
    }
    else
    {
    std::cout << "(";
    if(n.fanin[0].inv)
    {
      std::cout << "!";
    }
    if(n.type == node_type::mux_)
    {
      std::cout << "!" << char('a' + n.fanin[2].index) << " * ";
    }
    rec_print_dsd(rule, rule[n.fanin[0].index]);
    std::cout << " " << to_string(n.type) << " ";
    if(n.fanin[1].inv)
    {
      std::cout << "!";
    }
    if(n.type == node_type::mux_)
    {
      std::cout << char('a' + n.fanin[2].index) << " * ";
    }
    rec_print_dsd(rule, rule[n.fanin[1].index]);
    std::cout << ")";
    }
  }

  int try_top_dec(kitty::dynamic_truth_table& tt, int num_vars)
  {
    int i = 0;
    for(;i < num_vars;i++)
    {
      auto res = is_top_dec(tt, i);
      if(res.type != node_type::none)
        break;
    }
    return i;
  }

  dsd_node do_top_dec(kitty::dynamic_truth_table& tt, int index, std::vector<int> mapped_support)
  {
    auto node = is_top_dec(tt, index, &tt, false);

    node.fanin[0].index = mapped_support[index];
    return node;

  }

  std::tuple<int, int> try_bottom_dec(kitty::dynamic_truth_table& tt, int num_vars)
  {
    int i;
    int j;
    dsd_node res;
    for(i = 0;i < num_vars;i++)
    {
      for(j = i + 1; j < num_vars;j++)
      {
      res = is_bottom_dec(tt, i, j);
      if(res.type != node_type::none)
        break;
      }
      if(res.type != node_type::none)
        break;
    }
    std::tuple<int, int> ret = {i, j};
    return ret;
  }

  dsd_node do_bottom_dec(kitty::dynamic_truth_table& tt, int i, int j, int new_index, std::vector<int>& mapped_support)
  {
    auto node = is_bottom_dec(tt, i, j, &tt, new_index, false);

    node.fanin[0].index = mapped_support[i];
    node.fanin[1].index = mapped_support[j];

    mapped_support[i] = node.index;
    return node;
  }

  dsd_node do_shannon_dec(kitty::dynamic_truth_table tt, int index, kitty::dynamic_truth_table& co0, kitty::dynamic_truth_table& co1, std::vector<int> mapped_support)
  {
    auto node = shannon_dec(tt, index, &co0, &co1);
    node.fanin[0].index = mapped_support[index];
    return node;
  }

  void update_support(std::vector<int>& v, int index)
  {
    int i = 0;
    for(;i < v.size() && i < index; i++);

    for(;i < v.size(); i++)
    {
      v[i] = v[i+1];
    }
    v.pop_back();
  }

  template<class TT>
  void min_base_shrink(TT& tt, TT& tt_shr)
  {
    kitty::min_base_inplace(tt);
    kitty::shrink_to_inplace(tt_shr,tt);
  }

  int compute_rule(kitty::dynamic_truth_table& tt, std::vector<int> mapped_support, std::vector<dsd_node>& rule)
  {
    //tt has been already found
    if(!get_rules(tt).empty())
    {
      int count_old = 0;
      int count_curr = 0;
      std::vector<dsd_node> new_rule;
      auto found_rule = get_rules(tt);
      std::copy_if(found_rule.begin(), found_rule.end(), std::back_inserter(new_rule), [](dsd_node n) {
        return n.type != node_type::pi_;
      });
      for_each(found_rule.begin(), found_rule.end(), [&](dsd_node elem)
      {
        if(elem.type == node_type::pi_)
          count_old++;
      });
      for_each(rule.begin(), rule.end(), [&](dsd_node elem)
      {
          count_curr++;
      });
      //update index of node
      std::transform(new_rule.begin(), new_rule.end(), new_rule.begin(), [&](dsd_node& n) -> dsd_node
      {
        return {n.type, n.index + count_curr - count_old, n.fanin};
      });
      //update index of signal of fanins of nodes
      std::transform(new_rule.begin(), new_rule.end(), new_rule.begin(), [&](dsd_node& n) -> dsd_node
      {
        transform(n.fanin.begin(), n.fanin.end(), n.fanin.begin(), [&](signal s) -> signal
        {
          if(s.index >= count_old)
            return {s.inv, s.index + count_curr - count_old};
          else
            return {s.inv, mapped_support[s.index]};
        });       
        return {n.type, n.index, n.fanin};
      });
      rule.insert(rule.end(), new_rule.begin(), new_rule.end());
      return rule.size() - 1;
    }
    //top decomposition
    int i = try_top_dec(tt, tt.num_vars());
    if(i < tt.num_vars()) //it was top decomposable
    {
      auto res = do_top_dec(tt, i, mapped_support);

      update_support(mapped_support, i);

      kitty::dynamic_truth_table tt_shr(tt.num_vars()-1);
      min_base_shrink(tt, tt_shr);

      if(is_PI(tt_shr, tt_shr.num_vars()) < 0 && is_inv_PI(tt_shr, tt_shr.num_vars()) < 0) //check if remainder is PI
      {
        res.fanin.push_back({0,compute_rule(tt_shr, mapped_support, rule)});
      }
      else
      {
        if(is_PI(tt_shr, tt_shr.num_vars()) >= 0)
        {
          res.fanin.push_back({0, mapped_support[is_PI(tt_shr, tt_shr.num_vars())]});
        }
        else
        {
          res.fanin.push_back({1, mapped_support[is_inv_PI(tt_shr, tt_shr.num_vars())]});
        }
      }
      res.index = rule.size();      
      rule.push_back(res);
      return res.index;
    }
    else //bottom decomposition
    {
      auto couple = try_bottom_dec(tt, tt.num_vars());
      i = std::get<0>(couple);
      int j = std::get<1>(couple);

      if(i < tt.num_vars()) //it was bottom decomposable
      {
        auto res = do_bottom_dec(tt, i, j, rule.size(), mapped_support);
        rule.push_back(res);

        update_support(mapped_support, j);

        kitty::dynamic_truth_table tt_shr(tt.num_vars()-1);
        min_base_shrink(tt, tt_shr);

        return compute_rule(tt_shr, mapped_support, rule);
      }
      else //shannon dec
      {
        kitty::dynamic_truth_table co0(tt.num_vars());
        kitty::dynamic_truth_table co1(tt.num_vars());
        kitty::dynamic_truth_table co0_shr(tt.num_vars()-1);
        kitty::dynamic_truth_table co1_shr(tt.num_vars()-1);

        int index = find_unate_var(tt);

        auto res = do_shannon_dec(tt, index, co0, co1, mapped_support);
        
        //needed for PI checks
        int inv_var_co1 = is_inv_PI(co1, co1.num_vars());
        int map_inv_var_co1 = mapped_support[inv_var_co1];
        int var_co1 = is_PI(co1, co1.num_vars());
        int map_var_co1 = mapped_support[var_co1];
        int inv_var_co0 = is_inv_PI(co0, co0.num_vars());
        int map_inv_var_co0 = mapped_support[inv_var_co0];
        int var_co0 = is_PI(co0, co0.num_vars());
        int map_var_co0 = mapped_support[var_co0];

        update_support(mapped_support, index);

        if(inv_var_co1 < 0 && var_co1 < 0) //check if co1 is PI
        {
          min_base_shrink(co1, co1_shr);
          res.fanin.insert(res.fanin.begin(), {0, compute_rule(co1_shr, mapped_support, rule)});
        }
        else
        {
          if(inv_var_co1 >= 0)
          {
            res.fanin.insert(res.fanin.begin(), {1, map_inv_var_co1});
          }
          else
            res.fanin.insert(res.fanin.begin(), {0, map_var_co1});
        }
        if(inv_var_co0 < 0 && var_co0 < 0) //check if co0 is PI
        {
          min_base_shrink(co0, co0_shr);
          res.fanin.insert(res.fanin.begin(), {0, compute_rule(co0_shr, mapped_support, rule)});
        }
        else
        {
          if(inv_var_co0 >= 0)
          {
            res.fanin.insert(res.fanin.begin(), {1, map_inv_var_co0});
          }
          else
            res.fanin.insert(res.fanin.begin(), {0, map_var_co0});
        }      
        res.index = rule.size();      
        rule.push_back(res);
        return res.index;
      }
    }
  }

  int is_PI(kitty::dynamic_truth_table const& rem, int n_vars)
  {
    for(int i = 0; i < n_vars; i++)
    {
      auto var = rem.construct();
      kitty::create_nth_var( var, i );
      if(rem == var)
      {
        return i;
      }
    }
    return -1;
  }

  int is_inv_PI(kitty::dynamic_truth_table const& rem, int n_vars)
  {
    for(int i = 0; i < n_vars; i++)
    {
      auto var = rem.construct();
      kitty::create_nth_var( var, i );
      if(rem == ~var)
      {
        return i;
      }
    }
    return -1;
  }

private:

  template<class TT>  
  dsd_node is_top_dec( const TT& tt, int var_index, TT* func = nullptr, bool allow_xor = false )
  {
    static_assert( kitty::is_complete_truth_table<TT>::value, "Can only be applied on complete truth tables." );

    auto var = tt.construct();
    kitty::create_nth_var( var, var_index );

    if ( kitty::implies( tt, var ) )
    {
      if ( func )
      {
        *func = kitty::cofactor1( tt, var_index );
      }
      dsd_node res = {node_type::and_, var_index, {}};
      res.fanin.push_back({0, -1});
      return res;
    }
    else if ( kitty::implies( var, tt ) )
    {
      if ( func )
      {
        *func = kitty::cofactor0( tt, var_index );
      }
      dsd_node res = {node_type::or_, var_index, {}};
      res.fanin.push_back({0, -1});
      return res;
    }
    else if ( kitty::implies( tt, ~var ) )
    {
      if ( func )
      {
        *func = kitty::cofactor0( tt, var_index );
      }
      dsd_node res = {node_type::and_, var_index, {}};
      res.fanin.push_back({1, -1});
      return res;
    }
    else if ( kitty::implies( ~var, tt ) )
    {
      if ( func )
      {
        *func = kitty::cofactor1( tt, var_index );
      }
      dsd_node res = {node_type::or_, var_index, {}};
      res.fanin.push_back({1, -1});
      return res;
    }

    if ( allow_xor )
    {
      /* try XOR */
      const auto co0 = kitty::cofactor0( tt, var_index );
      const auto co1 = kitty::cofactor1( tt, var_index );

      if ( kitty::equal( co0, ~co1 ) )
      {
        if ( func )
        {
          *func = co0;
        }
        dsd_node res = {node_type::xor_, var_index, {}};
        res.fanin.push_back({0, -1});
        return res;
      }
    }

    return {node_type::none, var_index, {}};
  }

  template<class TT>
  dsd_node is_bottom_dec( const TT& tt, int var_index1, int var_index2, TT* func = nullptr, int new_index = -1, bool allow_xor = false )
  {
    static_assert( kitty::is_complete_truth_table<TT>::value, "Can only be applied on complete truth tables." );

    const auto tt0 = kitty::cofactor0( tt, var_index1 );
    const auto tt1 = kitty::cofactor1( tt, var_index1 );

    const auto tt00 = kitty::cofactor0( tt0, var_index2 );
    const auto tt01 = kitty::cofactor1( tt0, var_index2 );
    const auto tt10 = kitty::cofactor0( tt1, var_index2 );
    const auto tt11 = kitty::cofactor1( tt1, var_index2 );

    const auto eq01 = kitty::equal( tt00, tt01 );
    const auto eq02 = kitty::equal( tt00, tt10 );
    const auto eq03 = kitty::equal( tt00, tt11 );
    const auto eq12 = kitty::equal( tt01, tt10 );
    const auto eq13 = kitty::equal( tt01, tt11 );
    const auto eq23 = kitty::equal( tt10, tt11 );

    const auto num_pairs =
        static_cast<uint32_t>( eq01 ) +
        static_cast<uint32_t>( eq02 ) +
        static_cast<uint32_t>( eq03 ) +
        static_cast<uint32_t>( eq12 ) +
        static_cast<uint32_t>( eq13 ) +
        static_cast<uint32_t>( eq23 );

    if ( num_pairs != 2u && num_pairs != 3 )
    {
      return {node_type::none, -1, {}};
    }

    if ( !eq01 && !eq02 && !eq03 ) // 00 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt11, tt00 );
      }
      dsd_node res = {node_type::or_, new_index, {}};
      res.fanin.push_back({0, var_index1});
      res.fanin.push_back({0, var_index2});
      return res;
    }
    else if ( !eq01 && !eq12 && !eq13 ) // 01 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt10 );
      }
      dsd_node res = {node_type::and_, new_index, {}};
      res.fanin.push_back({1, var_index1});
      res.fanin.push_back({0, var_index2});
      return res;
    }
    else if ( !eq02 && !eq12 && !eq23 ) // 10 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt10 );
      }
      dsd_node res = {node_type::or_, new_index, {}};
      res.fanin.push_back({1, var_index1});
      res.fanin.push_back({0, var_index2});
      return res;
    }
    else if ( !eq03 && !eq13 && !eq23 ) // 11 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt11, tt00 );
      }
      dsd_node res = {node_type::and_, new_index, {}};
      res.fanin.push_back({0, var_index1});
      res.fanin.push_back({0, var_index2});
      return res;
    }
    else if ( allow_xor ) // XOR
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt00 );
      }
      dsd_node res = {node_type::xor_, new_index, {}};
      res.fanin.push_back({0, var_index1});
      res.fanin.push_back({0, var_index2});
      return res;
    }

    return {node_type::none, -1, {}};
  }

template<class TT>
int find_unate_var(const TT tt)
{
  int index;
    for(index = tt.num_vars() - 1; index > 0; index--)
    {
      const auto tt0 = kitty::cofactor0( tt, index );
      const auto tt1 = kitty::cofactor1( tt, index );
      if(((tt0 & tt1) == tt0) && ((tt0 & tt1) == tt1))
        return index;
    }
  return index;
}  

template<class TT>
dsd_node shannon_dec( const TT& tt, int index, TT* func0 = nullptr, TT* func1 = nullptr )
{
  static_assert( kitty::is_complete_truth_table<TT>::value, "Can only be applied on complete truth tables." );

    const auto tt0 = kitty::cofactor0( tt, index );
    const auto tt1 = kitty::cofactor1( tt, index );

    dsd_node res = {node_type::mux_, index, {}};
    res.fanin.push_back({0, index});

    if(func0 && func1)
    {
      *func0 = tt0;
      *func1 = tt1;
    }

    return res;
}

  bool _use_supergates; 

  std::vector<gate> const _gates; /* collection of gates */
  super_lib const& _supergates_spec; /* collection of supergates declarations */
  tech_library_params const _ps;
  super_utils<NInputs> _super; /* supergates generation */
  lib_t _super_lib; /* library of enumerated gates */
  lib_rule _rules_map; /*hash map for rules generation*/
};

}

/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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

/*!
  \file struct_library.hpp
  \brief Implements utilities for structural matching

  \author Alessandro Tempia Calvino
  \author Gianluca Radi
*/

#pragma once

#include <algorithm>
#include <cassert>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/decomposition.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include <parallel_hashmap/phmap.h>

#include "../io/genlib_reader.hpp"
#include "../io/super_reader.hpp"
#include "include/supergate.hpp"
#include "super_utils.hpp"
#include "tech_library.hpp"

namespace mockturtle
{

/*! \brief Library of gates for structural matching
 *
 * This class creates a technology library from a set
 * of input gates.
 *
 * Gates are processed to derive rules in the AIG format.
 * Then, every rule and subrule gets a unique id and the AND table is built.
 * Every gate gets a unique label comprehensive of its rule id and whether it is positive or negative.
 *
 * The template parameter `NInputs` selects the maximum number of variables
 * allowed for a gate in the library.
 *
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      std::vector<gate> gates;
      lorina::read_genlib( "file.genlib", genlib_reader( gates ) );
      // standard library
      mockturtle::struct_library lib( gates );
   \endverbatim
 */

template<unsigned NInputs = 9u>
class struct_library
{
public:
  enum class node_type
  {
    none,
    zero_,
    pi_,
    and_,
    or_,
    mux_,
    xor_
  };

  struct signal
  {
    union
    {
      struct
      {
        uint32_t inv : 1;
        uint32_t index : 31;
      };
      uint32_t data;
    };

    bool operator==( signal const& other ) const
    {
      return data == other.data;
    }
  };

  /*  struct for representing nodes in dsd decomposition */
  struct dsd_node
  {
    node_type type;

    int index;

    std::vector<signal> fanin = {};
  };

  /* struct for labels to assign to gates */
  struct label
  {
    union
    {
      struct
      {
        uint32_t inv : 1;
        uint32_t index : 31;
      };
      uint32_t data;
    };
    bool operator==( label const& other ) const
    {
      return data == other.data;
    }
  };

  struct signal_hash
  {
    std::size_t operator()( signal const& s ) const noexcept
    {
      return std::hash<uint32_t>{}( s.data );
    }
  };

  struct tuple_s_hash
  {
    std::size_t operator()( std::tuple<signal, signal> const& t ) const noexcept
    {
      size_t h1 = signal_hash()( std::get<0>( t ) );
      size_t h2 = signal_hash()( std::get<1>( t ) );
      return (uint64_t)h1 ^ ( ( (uint64_t)h2 ) << 32 ); // or use boost::hash_combine
    }
  };

private:
  using supergates_list_t = std::vector<supergate<NInputs>>;
  using composed_list_t = std::vector<composed_gate<NInputs>>;
  using lib_rule = phmap::flat_hash_map<kitty::dynamic_truth_table, std::vector<dsd_node>, kitty::hash<kitty::dynamic_truth_table>>;
  using rule = std::vector<dsd_node>;
  using lib_table = phmap::flat_hash_map<std::tuple<signal, signal>, uint32_t, tuple_s_hash>;
  using map_label_gate = std::unordered_map<uint32_t, supergates_list_t>;

public:
  explicit struct_library( std::vector<gate> const& gates )
      : _gates( gates ),
        _supergates(),
        _dsd_map(),
        _and_table(),
        _label_to_gate()
  {}

public:
  /*! \brief Construct the structural library.
   *
   * Generates the patterns for structural matching.
   * Variable `min_vars` defines the minimum number of
   * gate inputs considered for the library creation.
   * 0 < min_vars < UINT32_MAX
   */
  void construct( uint32_t min_vars = 2u, bool verbose = false, bool very_verbose = false )
  {
    generate_library( min_vars, verbose, very_verbose );
  }

  /*! \brief Construct the structural library.
   *
   * Generates the patterns for structural matching.
   */
  const map_label_gate& get_struct_library() const
  {
    return _label_to_gate;
  }

  /*! \brief Get the pattern ID.
   *
   *  \param id1 first pattern id.
   *  \param id2 second pattern id.
   * Returns a pattern ID if found, UINT32_MAX otherwise given the
   * children IDs. This function works with only AND operators.
   */
  const uint32_t get_pattern_id( uint32_t id1, uint32_t id2 ) const
  {
    signal l, r;
    l.data = id1;
    /* ignore input negations */
    if ( l.data == 3 )
      l.data = 2;
    r.data = id2;
    if ( r.data == 3 )
      r.data = 2;
    std::tuple<signal, signal> key;
    if ( l.index <= r.index )
      key = std::make_tuple( l, r );
    else
      key = std::make_tuple( r, l );
    auto match = _and_table.find( key );
    if ( match != _and_table.end() )
      return match->second;
    return UINT32_MAX;
  }

  /*! \brief Get the gates matching the pattern ID.
   *
   * Returns a list of gates that match the pattern ID.
   */
  const supergates_list_t* get_supergates_pattern( uint32_t id, bool phase ) const
  {
    auto match = _label_to_gate.find( ( id << 1 ) | ( phase ? 1 : 0 ) );
    if ( match != _label_to_gate.end() )
    {
      return &( match->second );
    }
    return nullptr;
  }

  /*! \brief Print and table.
   *
   */
  void print_and_table()
  {
    for ( auto elem : _and_table )
    {
      auto first0 = std::get<0>( elem.first );
      auto first1 = std::get<1>( elem.first );
      std::cout << "<" << ( first0.inv ? "!" : "" ) << first0.index;
      std::cout << ", " << ( first1.inv ? "!" : "" ) << first1.index << "> ";
      std::cout << elem.second << "\n";
    }
  }

private:
  void generate_library( uint32_t min_vars, bool verbose = false, bool very_verbose = false )
  {
    /* select and load gates */
    _supergates.reserve( _gates.size() );
    generate_composed_gates();

    /* mark dominate gates */
    std::vector<bool> skip_gates( _supergates.size(), false );
    // select_dominated_gates( skip_gates );

    std::vector<uint32_t> indexes( _supergates.size() );
    std::iota( indexes.begin(), indexes.end(), 0 );
    uint32_t max_label = 1;
    uint32_t gate_pol = 0; // polarity of AND equivalent gate
    uint32_t shift = 0;

    std::cout << "[i] processing structural library\n";
    if ( verbose )
      std::cout << "[i] structural library: processing " << _supergates.size() << " gates with minum size " << min_vars << " and max size " << NInputs << "\n";

    /* sort cells by increasing order of area */
    std::sort( indexes.begin(), indexes.end(),
               [&]( auto const& a, auto const& b ) -> bool {
                 return _supergates[a].area < _supergates[b].area;
               } );

    for ( uint32_t const ind : indexes )
    {
      composed_gate<NInputs> const& gate = _supergates[ind];

      if ( gate.num_vars < 2 || skip_gates[ind] )
        continue;

      /* DSD decomposition */
      rule rule = {};
      std::vector<int> support = {};
      for ( int i = 0; i < gate.num_vars; i++ )
      {
        rule.push_back( { node_type::pi_, i, {} } );
        support.push_back( i );
      }
      auto cpy = gate.function;
      gate_disjoint = false;
      compute_dsd( cpy, support, rule );

      if ( gate_disjoint )
        continue;

      _dsd_map.insert( { gate.function, rule } );
      if ( very_verbose )
      {
        std::cout << gate.root->name << "\n";
        std::cout << "Dsd:\n";
        print_rule( rule, rule[rule.size() - 1] );
      }

      /* Aig conversion */
      auto aig_rule = map_to_aig( rule );
      if ( very_verbose )
      {
        std::cout << "\nAig:\n";
        print_rule( aig_rule, aig_rule[aig_rule.size() - 1] );
      }

      /* Rules derivation */
      std::vector<std::vector<dsd_node>> der_rules = {};
      der_rules.push_back( aig_rule );
      std::vector<std::tuple<uint32_t, uint32_t>> depths = { { get_depth( aig_rule, aig_rule[aig_rule[aig_rule.size() - 1].fanin[0].index] ), get_depth( aig_rule, aig_rule[aig_rule[aig_rule.size() - 1].fanin[1].index] ) } };
      create_rules_from_dsd( der_rules, aig_rule, aig_rule[aig_rule.size() - 1], depths, true, true );

      /*std::vector<std::vector<dsd_node>> perm_rules = {};
      for( auto r : der_rules )
      {
        std::vector<uint32_t> variables = {};
        std::vector<uint32_t> variables_vis = {};
        std::vector<uint32_t> visited = {};
        uint32_t start_index = r.size()-1;
        std::vector<std::vector<uint32_t>> permut = {};

        while( variables_vis.size() < gate.num_vars )
        {
          variables.clear();

          find_batch( r, r[start_index], visited, variables );

          //debug
          std::cout << "visited\n";
          for(auto elem : visited)
            std::cout << elem << "\t";

          for(auto elem : variables)
            variables_vis.push_back(elem);

          if(variables.size() == 0)
          {
            start_index = find_start_index( r, r.size()-1, visited );
            continue;
          }

          permut = permute(variables);

          //first permutation already observed
          permut.erase(permut.begin());

           //debug
          std::cout << "perms\n";
          for(auto vec : permut)
          {
            for(auto elem : vec)
              std::cout << elem << "\t";
            std::cout << "\n";
          }

          if(permut.size() > 0)
            apply_perm(perm_rules, r, permut, variables);

          start_index = find_start_index( r, r.size()-1, visited );
        }
      }

      for( auto& rule : perm_rules )
      {
        der_rules.push_back(rule);
      }

      if ( very_verbose )
      {
        std::cout << "\nDerived:\n";
      }*/

      /* Indexing of rules and subrules, and_table construction, and gates' label assignement */
      for ( auto elem : der_rules )
      {
        gate_pol = 0;
        shift = 0;
        std::vector<uint8_t> perm( gate.num_vars );
        auto index_rule = do_indexing_rule( elem, elem[elem.size() - 1], max_label, gate_pol, perm, shift );

        /* skip gate creation for small gates (<`min_vars` inputs) and convergence */
        if ( gate.num_vars < min_vars )
        {
          continue;
          if ( very_verbose )
            std::cout << "[i] skip gate due to variable convergence\n";
        }

        // std::cout << mask_k << " " << mask_k1 << "\n";
        // std::cout << ((((gate_pol & ~(1 << k) ) & ~(1 << (k+1))) | (mask_k << 1)) | (mask_k1 >> 1)) << "\n";

        supergate<NInputs> sg = { &gate,
                                  static_cast<float>( gate.area ),
                                  gate.tdelay,
                                  perm,
                                  //((((gate_pol & ~(1 << k) ) & ~(1 << (k+1))) | (mask_k << 1)) | (mask_k1 >> 1)) };
                                  gate_pol };

        // std::cout << "Expression " << gate.root->expression << "\n";
        // std::cout << "polarity " << sg.polarity << "\n";

        for ( auto i = 0u; i < perm.size() && i < NInputs; ++i )
        {
          sg.tdelay[i] = gate.tdelay[sg.permutation[i]];
          // std::cout << "perm " << (uint32_t) sg.permutation[i] << "\t";
        }
        // std::cout << "\n";

        auto& v = _label_to_gate[index_rule.data];

        auto it = std::lower_bound( v.begin(), v.end(), sg, [&]( auto const& s1, auto const& s2 ) {
          if ( s1.area < s2.area )
            return true;
          if ( s1.area > s2.area )
            return false;
          if ( s1.root->num_vars < s2.root->num_vars )
            return true;
          if ( s1.root->num_vars > s2.root->num_vars )
            return true;
          return s1.root->id < s2.root->id;
        } );

        v.insert( it, sg );

        if ( very_verbose )
        {
          print_rule( elem, elem[elem.size() - 1] );
          std::cout << "\n\n";
        }
      }
    }
    if ( very_verbose )
      std::cout << "\n";
  }

  void generate_composed_gates()
  {
    /* filter multi-output gates */
    std::unordered_map<std::string, uint32_t> multioutput_map;
    multioutput_map.reserve( _gates.size() );

    for ( const auto& g : _gates )
    {
      if ( multioutput_map.find( g.name ) != multioutput_map.end() )
      {
        multioutput_map[g.name] += 1;
      }
      else
      {
        multioutput_map[g.name] = 1;
      }
    }

    /* create composed gates */
    uint32_t ignored = 0;
    for ( const auto& g : _gates )
    {
      std::array<float, NInputs> pin_to_pin_delays{};

      /* filter large gates and multi-output gates */
      if ( g.function.num_vars() > NInputs || multioutput_map[g.name] > 1 )
      {
        ++ignored;
        continue;
      }

      auto i = 0u;
      for ( auto const& pin : g.pins )
      {
        /* use worst pin delay */
        pin_to_pin_delays[i++] = std::max( pin.rise_block_delay, pin.fall_block_delay );
      }

      _supergates.emplace_back( composed_gate<NInputs>{ static_cast<unsigned int>( _supergates.size() ),
                                                        false,
                                                        &g,
                                                        g.num_vars,
                                                        g.function,
                                                        g.area,
                                                        pin_to_pin_delays,
                                                        {} } );
    }
  }

  void select_dominated_gates( std::vector<bool>& skip_gates )
  {
    for ( uint32_t i = 0; i < skip_gates.size() - 1; ++i )
    {
      if ( _supergates[i].root == nullptr )
        continue;

      if ( skip_gates[i] )
        continue;

      auto const& tti = _supergates[i].function;
      for ( uint32_t j = i + 1; j < skip_gates.size(); ++j )
      {
        auto const& ttj = _supergates[j].function;

        /* get the same functionality */
        if ( tti != ttj )
          continue;

        /* is i smaller than j */
        bool smaller = _supergates[i].area < _supergates[j].area;

        /* is i faster for every pin */
        bool faster = true;
        for ( uint32_t k = 0; k < tti.num_vars(); ++k )
        {
          if ( _supergates[i].tdelay[k] > _supergates[j].tdelay[k] )
            faster = false;
        }

        if ( smaller && faster )
        {
          skip_gates[j] = true;
          continue;
        }

        /* is j faster for every pin */
        faster = true;
        for ( uint32_t k = 0; k < tti.num_vars(); ++k )
        {
          if ( _supergates[j].tdelay[k] > _supergates[i].tdelay[k] )
            faster = false;
        }

        if ( !smaller && faster )
        {
          skip_gates[i] = true;
          break;
        }
      }
    }
  }

  int try_top_dec( kitty::dynamic_truth_table& tt, int num_vars )
  {
    int i = 0;
    for ( ; i < num_vars; i++ )
    {
      auto res = is_top_dec( tt, i, false );
      if ( res.type != node_type::none )
        break;
    }
    return i;
  }

  dsd_node do_top_dec( kitty::dynamic_truth_table& tt, int index, std::vector<int> mapped_support )
  {
    auto node = is_top_dec( tt, index, false, &tt );

    node.fanin[0].index = mapped_support[index];
    return node;
  }

  std::tuple<int, int> try_bottom_dec( kitty::dynamic_truth_table& tt, int num_vars )
  {
    int i;
    int j;
    dsd_node res;
    for ( i = 0; i < num_vars; i++ )
    {
      for ( j = i + 1; j < num_vars; j++ )
      {
        res = is_bottom_dec( tt, i, j );
        if ( res.type != node_type::none )
          break;
      }
      if ( res.type != node_type::none )
        break;
    }
    std::tuple<int, int> ret = { i, j };
    return ret;
  }

  dsd_node do_bottom_dec( kitty::dynamic_truth_table& tt, int i, int j, int new_index, std::vector<int>& mapped_support )
  {
    auto node = is_bottom_dec( tt, i, j, &tt, new_index, false );

    node.fanin[0].index = mapped_support[i];
    node.fanin[1].index = mapped_support[j];

    mapped_support[i] = node.index;
    return node;
  }

  dsd_node do_shannon_dec( kitty::dynamic_truth_table tt, int index, kitty::dynamic_truth_table& co0, kitty::dynamic_truth_table& co1, std::vector<int> mapped_support )
  {
    auto node = shannon_dec( tt, index, &co0, &co1 );
    node.fanin[0].index = mapped_support[index];
    return node;
  }

  void update_support( std::vector<int>& v, int index )
  {
    int i = 0;
    for ( ; i < v.size() && i < index; i++ )
      ;

    for ( ; i < v.size(); i++ )
    {
      v[i] = v[i + 1];
    }
    v.pop_back();
  }

  template<class TT>
  void min_base_shrink( TT& tt, TT& tt_shr )
  {
    kitty::min_base_inplace( tt );
    kitty::shrink_to_inplace( tt_shr, tt );
  }

  int is_PI( kitty::dynamic_truth_table const& rem, int n_vars )
  {
    for ( int i = 0; i < n_vars; i++ )
    {
      auto var = rem.construct();
      kitty::create_nth_var( var, i );
      if ( rem == var )
      {
        return i;
      }
    }
    return -1;
  }

  int is_inv_PI( kitty::dynamic_truth_table const& rem, int n_vars )
  {
    for ( int i = 0; i < n_vars; i++ )
    {
      auto var = rem.construct();
      kitty::create_nth_var( var, i );
      if ( rem == ~var )
      {
        return i;
      }
    }
    return -1;
  }

  void update_found_rule( kitty::dynamic_truth_table& tt, std::vector<int>& mapped_support, std::vector<dsd_node>& rule )
  {
    int count_old = 0;
    int count_curr = 0;
    std::vector<dsd_node> new_rule;
    auto found_rule = get_rules( tt );
    std::copy_if( found_rule.begin(), found_rule.end(), std::back_inserter( new_rule ), []( dsd_node n ) {
      return ( n.type != node_type::pi_ );
    } );
    for_each( found_rule.begin(), found_rule.end(), [&]( dsd_node elem ) {
      if ( elem.type == node_type::pi_ )
        count_old++;
    } );
    for_each( rule.begin(), rule.end(), [&]( dsd_node elem ) {
      count_curr++;
    } );
    /* update index of node */
    std::transform( new_rule.begin(), new_rule.end(), new_rule.begin(), [&]( dsd_node& n ) -> dsd_node {
      return { n.type, n.index + count_curr - count_old, n.fanin };
    } );
    /* update index of signal of fanins of nodes */
    std::transform( new_rule.begin(), new_rule.end(), new_rule.begin(), [&]( dsd_node& n ) -> dsd_node {
      transform( n.fanin.begin(), n.fanin.end(), n.fanin.begin(), [&]( signal s ) -> signal {
        if ( s.index >= count_old )
          return { s.inv, s.index + count_curr - count_old };
        else
          return { s.inv, mapped_support[s.index] };
      } );
      return { n.type, n.index, n.fanin };
    } );
    rule.insert( rule.end(), new_rule.begin(), new_rule.end() );
  }

  /*! \brief Compute DSD decomposition for a boolean function recursively.
   *
   *  \param tt dynamic truth table representing the function.
   *  \param mapped_support vector indicating function's support at every recursive step.
   *  \param rule DSD decomposition of the function.
   * Returns index of dsd_node to add to rule.
   */
  int compute_dsd( kitty::dynamic_truth_table& tt, std::vector<int> mapped_support, std::vector<dsd_node>& rule )
  {
    /* Function has been already found */
    if ( !get_rules( tt ).empty() )
    {
      update_found_rule( tt, mapped_support, rule );
      return rule.size() - 1;
    }
    /* try top decomposition */
    int i = try_top_dec( tt, tt.num_vars() );
    if ( i < tt.num_vars() ) // it was top decomposable
    {
      auto res = do_top_dec( tt, i, mapped_support );

      update_support( mapped_support, i );

      kitty::dynamic_truth_table tt_shr( tt.num_vars() - 1 );
      min_base_shrink( tt, tt_shr );

      if ( is_PI( tt_shr, tt_shr.num_vars() ) < 0 && is_inv_PI( tt_shr, tt_shr.num_vars() ) < 0 ) // check if remainder is PI
      {
        res.fanin.push_back( { 0, compute_dsd( tt_shr, mapped_support, rule ) } );
      }
      else
      {
        if ( is_PI( tt_shr, tt_shr.num_vars() ) >= 0 )
        {
          res.fanin.push_back( { 0, mapped_support[is_PI( tt_shr, tt_shr.num_vars() )] } );
        }
        else
        {
          res.fanin.push_back( { 1, mapped_support[is_inv_PI( tt_shr, tt_shr.num_vars() )] } );
        }
      }
      res.index = rule.size();
      rule.push_back( res );
      return res.index;
    }

    /* try bottom decomposition */
    else
    {
      auto couple = try_bottom_dec( tt, tt.num_vars() );
      i = std::get<0>( couple );
      int j = std::get<1>( couple );

      if ( i < tt.num_vars() ) // it was bottom decomposable
      {
        auto res = do_bottom_dec( tt, i, j, rule.size(), mapped_support );
        rule.push_back( res );

        update_support( mapped_support, j );

        kitty::dynamic_truth_table tt_shr( tt.num_vars() - 1 );
        min_base_shrink( tt, tt_shr );

        return compute_dsd( tt_shr, mapped_support, rule );
      }

      /* do shannon decomposition */
      else
      {
        kitty::dynamic_truth_table co0( tt.num_vars() );
        kitty::dynamic_truth_table co1( tt.num_vars() );
        kitty::dynamic_truth_table co0_shr( tt.num_vars() - 1 );
        kitty::dynamic_truth_table co1_shr( tt.num_vars() - 1 );

        int index = find_unate_var( tt );

        auto res = do_shannon_dec( tt, index, co0, co1, mapped_support );

        /* check for reconvergence */
        gate_disjoint = true;

        int inv_var_co1 = is_inv_PI( co1, co1.num_vars() );
        int map_inv_var_co1 = mapped_support[inv_var_co1];
        int var_co1 = is_PI( co1, co1.num_vars() );
        int map_var_co1 = mapped_support[var_co1];
        int inv_var_co0 = is_inv_PI( co0, co0.num_vars() );
        int map_inv_var_co0 = mapped_support[inv_var_co0];
        int var_co0 = is_PI( co0, co0.num_vars() );
        int map_var_co0 = mapped_support[var_co0];

        update_support( mapped_support, index );

        if ( inv_var_co1 < 0 && var_co1 < 0 ) // check if co1 is PI
        {
          min_base_shrink( co1, co1_shr );
          res.fanin.insert( res.fanin.begin(), { 0, compute_dsd( co1_shr, mapped_support, rule ) } );
        }
        else
        {
          if ( inv_var_co1 >= 0 )
          {
            res.fanin.insert( res.fanin.begin(), { 1, map_inv_var_co1 } );
          }
          else
            res.fanin.insert( res.fanin.begin(), { 0, map_var_co1 } );
        }
        if ( inv_var_co0 < 0 && var_co0 < 0 ) // check if co0 is PI
        {
          min_base_shrink( co0, co0_shr );
          res.fanin.insert( res.fanin.begin(), { 0, compute_dsd( co0_shr, mapped_support, rule ) } );
        }
        else
        {
          if ( inv_var_co0 >= 0 )
          {
            res.fanin.insert( res.fanin.begin(), { 1, map_inv_var_co0 } );
          }
          else
            res.fanin.insert( res.fanin.begin(), { 0, map_var_co0 } );
        }
        res.index = rule.size();
        rule.push_back( res );
        return res.index;
      }
    }
  }

  rule get_rules( kitty::dynamic_truth_table const& tt )
  {
    auto match = _dsd_map.find( tt );
    if ( match != _dsd_map.end() )
      return match->second;
    return {};
  }

  dsd_node* get_father( rule& rule, dsd_node& node )
  {
    for ( uint32_t i = 0; i < rule.size(); i++ )
    {
      if ( rule[i].type != node_type::pi_ && rule[i].type != node_type::zero_ && ( rule[i].fanin[0].index == node.index || rule[i].fanin[1].index == node.index ) )
        return &rule[i];
    }
    return NULL;
  }

  dsd_node* find_node( rule& r, uint32_t i )
  {
    for ( int j = 0; j < r.size(); j++ )
    {
      if ( r[j].index == i )
        return &r[j];
    }
    return NULL;
  }

  /*! \brief Convert rule derived from DSD decomposition into aig format.
   *
   *  \param r rule to convert.
   * Returns rule converted into aig format.
   */
  rule map_to_aig( rule& r )
  {
    std::vector<dsd_node> rule( r );
    std::vector<dsd_node> aig_rule;

    std::transform( rule.begin(), rule.end(), rule.begin(), []( dsd_node n ) -> dsd_node {
      for ( auto& s : n.fanin )
      {
        s.index += 1;
      }
      return { n.type, n.index + 1, n.fanin };
    } );

    rule.insert( rule.begin(), { node_type::zero_, 0, {} } );

    for ( typename std::vector<dsd_node>::reverse_iterator i = rule.rbegin(); i != rule.rend(); ++i )
    {
      dsd_node n = *i;
      dsd_node new_node;

      if ( n.type == node_type::and_ || n.type == node_type::pi_ || n.type == node_type::zero_ )
      {
        new_node = n;
      }
      else if ( n.type == node_type::or_ )
      {
        new_node = { node_type::and_, n.index, { { ~n.fanin[0].inv, n.fanin[0].index }, { ~n.fanin[1].inv, n.fanin[1].index } } };
        if ( get_father( rule, n ) != NULL )
        {
          dsd_node* father = find_node( aig_rule, get_father( rule, n )->index );
          if ( father->fanin[0].index == n.index )
          {
            father->fanin[0].inv = ~father->fanin[0].inv;
          }
          else
          {
            father->fanin[1].inv = ~father->fanin[1].inv;
          }
        }
        else // it is root
        {
          dsd_node new_root = { node_type::and_, rule.size(), { { 1, 0 }, { 1, n.index } } };
          aig_rule.insert( aig_rule.begin(), new_root );
        }
      }
      else if ( n.type == node_type::mux_ )
      {
        if ( get_father( rule, n ) != NULL )
        {
          dsd_node* father = find_node( aig_rule, get_father( rule, n )->index );
          if ( father->fanin[0].index == n.index )
          {
            father->fanin[0].inv = ~father->fanin[0].inv;
          }
          else
          {
            father->fanin[1].inv = ~father->fanin[1].inv;
          }
          dsd_node node_or = { node_type::and_, n.index + 2, { { 1, n.index }, { 1, n.index + 1 } } };
          dsd_node node_and1 = { node_type::and_, n.index + 1, { { 0, n.fanin[2].index }, { n.fanin[1].inv, n.fanin[1].index } } };
          new_node = { node_type::and_, n.index, { { 1, n.fanin[2].index }, { n.fanin[0].inv, n.fanin[0].index } } }; // and0_node

          // node already in aig_rule must have index and fanin index update (index += 2, fanin_index -> (>= n.index -> +2; nothing))
          for ( auto& elem : aig_rule )
          {
            elem.index += 2;
            for ( auto& s : elem.fanin )
            {
              if ( s.index >= n.index )
                s.index += 2;
            }
          }

          aig_rule.insert( aig_rule.begin(), node_or );
          aig_rule.insert( aig_rule.begin(), node_and1 );
        }
        else // it is root
        {
          dsd_node new_root = { node_type::and_, rule.size() + 2, { { 1, 0 }, { 1, rule.size() + 1 } } };
          dsd_node node_or = { node_type::and_, rule.size() + 1, { { 1, rule.size() - 1 }, { 1, rule.size() } } };
          dsd_node node_and1 = { node_type::and_, rule.size(), { { 0, n.fanin[2].index }, { n.fanin[1].inv, n.fanin[1].index } } };
          new_node = { node_type::and_, rule.size() - 1, { { 1, n.fanin[2].index }, { n.fanin[0].inv, n.fanin[0].index } } }; // and0_node

          aig_rule.insert( aig_rule.begin(), new_root );
          aig_rule.insert( aig_rule.begin(), node_or );
          aig_rule.insert( aig_rule.begin(), node_and1 );
        }
      }
      else if ( n.type == node_type::xor_ )
      {
        if ( get_father( rule, n ) != NULL )
        {
          dsd_node* father = find_node( aig_rule, get_father( rule, n )->index );
          if ( father->fanin[0].index == n.index )
          {
            father->fanin[0].inv = ~father->fanin[0].inv;
          }
          else
          {
            father->fanin[1].inv = ~father->fanin[1].inv;
          }
          dsd_node node_or = { node_type::and_, n.index + 2, { { 1, n.index }, { 1, n.index + 1 } } };
          dsd_node node_and1 = { node_type::and_, n.index + 1, { { 0, n.fanin[0].index }, { 1, n.fanin[1].index } } };
          new_node = { node_type::and_, n.index, { { 1, n.fanin[0].index }, { 0, n.fanin[1].index } } }; // and0_node

          // node already in aig_rule must have index and fanin index update (index += 2, fanin_index -> (>= n.index -> +2; nothing))
          for ( auto& elem : aig_rule )
          {
            elem.index += 2;
            for ( auto& s : elem.fanin )
            {
              if ( s.index >= n.index )
                s.index += 2;
            }
          }

          aig_rule.insert( aig_rule.begin(), node_or );
          aig_rule.insert( aig_rule.begin(), node_and1 );
        }
        else // it is root
        {
          dsd_node new_root = { node_type::and_, rule.size() + 2, { { 1, 0 }, { 1, rule.size() + 1 } } };
          dsd_node node_or = { node_type::and_, rule.size() + 1, { { 1, rule.size() - 1 }, { 1, rule.size() } } };
          dsd_node node_and1 = { node_type::and_, rule.size(), { { 0, n.fanin[0].index }, { 1, n.fanin[1].index } } };
          new_node = { node_type::and_, rule.size() - 1, { { 1, n.fanin[0].index }, { 0, n.fanin[1].index } } }; // and0_node

          aig_rule.insert( aig_rule.begin(), new_root );
          aig_rule.insert( aig_rule.begin(), node_or );
          aig_rule.insert( aig_rule.begin(), node_and1 );
        }
      }
      aig_rule.insert( aig_rule.begin(), new_node );
    }
    return aig_rule;
  }

  void swap( rule& rule, dsd_node* node_i, dsd_node* node_j )
  {
    auto i = node_i->index;
    auto j = node_j->index;
    node_i->index = j;
    node_j->index = i;
    std::swap( rule[i], rule[j] );
  }

  /* makes left or right move to derive a new rule */
  void make_move( rule& rule, dsd_node* target, dsd_node* r, uint8_t left )
  {
    auto targ_index = 0 + left;
    auto r_index = 1 - targ_index;
    auto temp_index = target->fanin[targ_index].index;
    auto temp_inv = target->fanin[targ_index].inv;
    // swap position internal to rule of the two elements
    swap( rule, target, r );
    auto temp = target;
    target = r;
    r = temp;

    // adjust children
    target->fanin[targ_index].index = r->index;
    target->fanin[targ_index].inv = 0;
    r->fanin[r_index].index = temp_index;
    r->fanin[r_index].inv = temp_inv;
  }

  /* checks whether a rule with the same left depth or right depth has already been encountered */
  bool check_depths( rule rule, dsd_node root, std::vector<std::tuple<uint32_t, uint32_t>> depth_branches, uint32_t left ) // for left = 1 check if left move is possible
  {
    auto left_node = rule[root.fanin[0].index];
    auto right_node = rule[root.fanin[1].index];
    auto left_depth = get_depth( rule, left_node );
    auto right_depth = get_depth( rule, right_node );
    auto left_it1 = std::find( depth_branches.begin(), depth_branches.end(), ( std::tuple<uint32_t, uint32_t> ){ left_depth - 1, right_depth + 1 } );
    auto left_it2 = std::find( depth_branches.begin(), depth_branches.end(), ( std::tuple<uint32_t, uint32_t> ){ right_depth + 1, left_depth - 1 } );
    auto right_it1 = std::find( depth_branches.begin(), depth_branches.end(), ( std::tuple<uint32_t, uint32_t> ){ left_depth + 1, right_depth - 1 } );
    auto right_it2 = std::find( depth_branches.begin(), depth_branches.end(), ( std::tuple<uint32_t, uint32_t> ){ right_depth - 1, left_depth + 1 } );
    if ( left )
      return left_it1 == depth_branches.end() && left_it2 == depth_branches.end();
    else
      return right_it1 == depth_branches.end() && right_it2 == depth_branches.end();
  }

  /*! \brief Recursively create new rules from the given one.
   *         For every dsd node of the original rule, left and right moves are tried.
   *         Then, the same algorithm is applied to all derived rules.
   *         The algorithm stops if no new acceptable rules can be found.
   *         A new rule is acceptable if no other rule with the same left and right depths has been found.
   *  \param new_rules vector of derived rules.
   *  \param rule original rule.
   *  \param start_node node of rule on which we try the right and left moves.
   *  \param depth_branches vector of encountered left and right depths.
   *  \param can_left specifies if we can perform a left move.
   *  \param can_right specifies if we can perform a right move.
   */
  void create_rules_from_dsd( std::vector<rule>& new_rules, rule rule, dsd_node start_node, std::vector<std::tuple<uint32_t, uint32_t>>& depth_branches, bool can_left, bool can_right )
  {
    if ( start_node.type == node_type::pi_ || start_node.type == node_type::zero_ ) // if you cannot produce new rules or you are a PI return
      return;

    std::vector<dsd_node> left_rule( rule );
    std::vector<dsd_node> right_rule( rule );
    std::vector<std::tuple<uint32_t, uint32_t>> next_depths = {};
    auto left_node = &left_rule[start_node.fanin[0].index];
    auto right_node = &right_rule[start_node.fanin[1].index];
    bool new_left = false;
    bool new_right = false;

    std::tuple<uint32_t, uint32_t> depths = { get_depth( rule, *left_node ), get_depth( rule, *right_node ) };
    depth_branches.push_back( depths );

    /* left move */
    if ( can_left && left_node->type == start_node.type && start_node.fanin[0].inv == 0 && ( !not_simmetry || check_depths( rule, start_node, depth_branches, 1 ) ) )
    {
      auto r = &left_rule[start_node.index];

      make_move( left_rule, left_node, r, 1 );

      new_rules.push_back( left_rule );
      new_left = true;
      depth_branches.push_back( { get_depth( left_rule, left_rule[left_rule[left_node->index].fanin[0].index] ), get_depth( left_rule, left_rule[left_rule[left_node->index].fanin[1].index] ) } );
    }
    /* right move */
    if ( can_right && right_node->type == start_node.type && start_node.fanin[1].inv == 0 && ( !not_simmetry || check_depths( rule, start_node, depth_branches, 0 ) ) )
    {
      auto r = &right_rule[start_node.index];

      make_move( right_rule, right_node, r, 0 );

      new_rules.push_back( right_rule );
      new_right = true;
      depth_branches.push_back( { get_depth( right_rule, right_rule[right_rule[right_node->index].fanin[0].index] ), get_depth( right_rule, right_rule[right_rule[right_node->index].fanin[1].index] ) } );
    }
    if ( permutation && left_node->type == node_type::pi_ && right_node->type == node_type::pi_ )
    {
      auto r = &right_rule[start_node.index];
      r->fanin[0].index = right_node->index;
      r->fanin[1].index = left_node->index;
      new_rules.push_back( right_rule );
    }

    /* initial rule, start_node left children */
    create_rules_from_dsd( new_rules, rule, rule[start_node.fanin[0].index], next_depths, true, true );
    /* initial rule, start_node right children */
    create_rules_from_dsd( new_rules, rule, rule[start_node.fanin[1].index], next_depths, true, true );
    /* left rule, start_node new root */
    if ( new_left )
    {
      create_rules_from_dsd( new_rules, left_rule, left_rule[start_node.index], depth_branches, true, false );
    }
    /* right rule, start_node new root */
    if ( new_right )
    {
      create_rules_from_dsd( new_rules, right_rule, right_rule[start_node.index], depth_branches, false, true );
    }
  }

  void permutations( std::vector<std::vector<uint32_t>>& res, std::vector<uint32_t> nums, uint32_t l, uint32_t h )
  {
    // Base case
    // Add the vector to result and return
    if ( l == h )
    {
      res.push_back( nums );
      return;
    }

    // Permutations made
    for ( uint32_t i = l; i <= h; i++ )
    {

      // Swapping
      std::swap( nums[l], nums[i] );

      // Calling permutations for
      // next greater value of l
      permutations( res, nums, l + 1, h );

      // Backtracking
      std::swap( nums[l], nums[i] );
    }
  }

  // Function to get the permutations
  std::vector<std::vector<uint32_t>> permute( std::vector<uint32_t>& nums )
  {
    // Declaring result variable
    std::vector<std::vector<uint32_t>> res;
    uint32_t x = nums.size() - 1;

    // Calling permutations for the first
    // time by passing l
    // as 0 and h = nums.size()-1
    permutations( res, nums, 0, x );
    return res;
  }

  void find_batch( rule rule, dsd_node start, std::vector<uint32_t>& visited, std::vector<uint32_t>& variables )
  {
    visited.push_back( start.index );
    if ( rule[start.fanin[0].index].type == node_type::pi_ )
    {
      variables.push_back( start.fanin[0].index );
    }
    if ( rule[start.fanin[1].index].type == node_type::pi_ )
    {
      variables.push_back( start.fanin[1].index );
    }
    if ( start.fanin[0].inv == 0 && rule[start.fanin[0].index].type != node_type::pi_ && std::find( visited.begin(), visited.end(), start.fanin[0].index ) == visited.end() )
    {
      find_batch( rule, rule[start.fanin[0].index], visited, variables );
    }
    if ( start.fanin[1].inv == 0 && rule[start.fanin[1].index].type != node_type::pi_ && std::find( visited.begin(), visited.end(), start.fanin[1].index ) == visited.end() )
    {
      find_batch( rule, rule[start.fanin[1].index], visited, variables );
    }
    return;
  }

  uint32_t find_start_index( rule r, uint32_t start_index, std::vector<uint32_t> visited )
  {
    auto start_node = r[start_index];
    if ( std::find( visited.begin(), visited.end(), start_index ) == visited.end() )
      return start_index;
    if ( r[start_node.fanin[0].index].type != node_type::pi_ && r[start_node.fanin[0].index].type != node_type::zero_ )
    {
      uint32_t res = find_start_index( r, r[start_index].fanin[0].index, visited );
      if ( res != std::numeric_limits<uint32_t>::max() )
        return res;
    }
    if ( r[start_node.fanin[1].index].type != node_type::pi_ && r[start_node.fanin[1].index].type != node_type::zero_ )
    {
      uint32_t res = find_start_index( r, r[start_index].fanin[1].index, visited );
      if ( res != std::numeric_limits<uint32_t>::max() )
        return res;
    }
    return std::numeric_limits<uint32_t>::max();
  }

  void apply_perm( std::vector<rule>& perm_rules, rule r, std::vector<std::vector<uint32_t>> perms, std::vector<uint32_t> variables )
  {
    for ( auto perm : perms )
    {
      std::vector<dsd_node> new_rule( r );
      for ( int i = 0; i < perm.size(); i++ )
      {
        dsd_node* father_new = get_father( new_rule, new_rule[variables[i]] );
        dsd_node* father = get_father( r, r[variables[i]] );
        if ( father->fanin[0].index == variables[i] )
        {
          father_new->fanin[0].index = perm[i];
          std::cout << "father left " << father_new->fanin[0].index << "\n";
        }
        else
        {
          father_new->fanin[1].index = perm[i];
          std::cout << "father right " << father_new->fanin[1].index << "\n";
        }
      }
      perm_rules.push_back( new_rule );
    }
  }

  uint32_t compute_canonized_polarity( uint32_t polarity, int left_pi, int right_pi, int obs_pi )
  {
    uint32_t mask_l = 0;
    uint32_t mask_r = 0;
    uint32_t mask_obs = 0;
    for ( int i = 0; i < obs_pi; i++ )
    {
      mask_obs |= ( 1 << i );
    }
    for ( int i = obs_pi; i < obs_pi + left_pi; i++ )
    {
      mask_l |= ( 1 << i );
    }
    for ( int i = obs_pi + left_pi; i < obs_pi + left_pi + right_pi; i++ )
    {
      mask_r |= ( 1 << i );
    }
    return ( ( polarity & mask_l ) << right_pi ) | ( ( polarity & mask_r ) >> left_pi ) | ( polarity & mask_obs );
  }

  void compute_canonized_permutation( std::vector<uint8_t>& perm, int left_pi, int right_pi, int obs_pi )
  {
    std::vector copy( perm );
    for ( int i = obs_pi; i < obs_pi + right_pi; i++ )
    {
      perm[i] = copy[( i + left_pi )];
    }
    for ( int i = right_pi + obs_pi; i < perm.size(); i++ )
    {
      perm[i] = copy[( i - right_pi )];
    }
  }

  /*! \brief Recursively assigns indexes to a rule and its subrules and builds and_table.
   *         It also computes negations and permutations for the gate whose rule is being passed as parameter.
   *
   *  \param r rule to index.
   *  \param n dsd node to start from.
   *  \param max max index assigned.
   *  \param polarity polarity of gate.
   *  \param perm permutation of gate.
   *  \param shift specifies the number of PIs encountered.
   * Returns label to be assigned to gate whose rule is r.
   */
  label do_indexing_rule( rule r, dsd_node n, uint32_t& max, uint32_t& polarity, std::vector<uint8_t>& perm, uint32_t& shift )
  {
    if ( n.type == node_type::pi_ )
    {
      perm[shift] = n.index - 1;
      return { 0, 1 };
    }
    if ( n.type == node_type::zero_ )
      return { 0, 0 };

    int obs_pi = shift;

    /* do indexing on the left */
    uint32_t left_index = do_indexing_rule( r, r[n.fanin[0].index], max, polarity, perm, shift ).index;
    if ( r[n.fanin[0].index].type == node_type::pi_ )
    {
      polarity |= ( n.fanin[0].inv << shift );
      shift++;
    }
    /* encountered PIs on the left */
    int left_pi = shift - obs_pi;

    /* do indexing on the right */
    uint32_t right_index = do_indexing_rule( r, r[n.fanin[1].index], max, polarity, perm, shift ).index;
    if ( r[n.fanin[1].index].type == node_type::pi_ )
    {
      polarity |= ( n.fanin[1].inv << shift );
      shift++;
    }
    /* encountered PIs on the right */
    int right_pi = shift - obs_pi - left_pi;

    /* check if it is inverted gate */
    if ( n.fanin[0].index == 0 && n.fanin[1].inv && n.index == r.size() - 1 )
    {
      return { 1, right_index };
    }
    signal left, right;

    /* ignore invertion of PIs */
    if ( r[n.fanin[0].index].type == node_type::pi_ )
      left.inv = 0;
    else
      left.inv = (uint64_t)n.fanin[0].inv;
    left.index = left_index;
    if ( r[n.fanin[1].index].type == node_type::pi_ )
      right.inv = 0;
    else
      right.inv = (uint64_t)n.fanin[1].inv;
    right.index = right_index;

    std::tuple<signal, signal> t;

    /* canonize and_table on left index being smaller than right one */
    if ( left.index <= right.index )
      t = std::make_tuple( left, right );
    else
    {
      /* new polarity */
      polarity = compute_canonized_polarity( polarity, left_pi, right_pi, obs_pi );

      /* new permutation */
      compute_canonized_permutation( perm, left_pi, right_pi, obs_pi );

      t = std::make_tuple( right, left );
    }
    auto match = _and_table.find( t );
    if ( match != _and_table.end() )
      return { 0, match->second };
    max++;
    /* insert new value in and_table */
    _and_table.insert( { t, max } );
    return { 0, max };
  }

  template<class TT>
  dsd_node is_top_dec( const TT& tt, int var_index, bool allow_xor = false, TT* func = nullptr )
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
      dsd_node res = { node_type::and_, var_index, {} };
      res.fanin.push_back( { 0, -1 } );
      return res;
    }
    else if ( kitty::implies( var, tt ) )
    {
      if ( func )
      {
        *func = kitty::cofactor0( tt, var_index );
      }
      dsd_node res = { node_type::or_, var_index, {} };
      res.fanin.push_back( { 0, -1 } );
      return res;
    }
    else if ( kitty::implies( tt, ~var ) )
    {
      if ( func )
      {
        *func = kitty::cofactor0( tt, var_index );
      }
      dsd_node res = { node_type::and_, var_index, {} };
      res.fanin.push_back( { 1, -1 } );
      return res;
    }
    else if ( kitty::implies( ~var, tt ) )
    {
      if ( func )
      {
        *func = kitty::cofactor1( tt, var_index );
      }
      dsd_node res = { node_type::or_, var_index, {} };
      res.fanin.push_back( { 1, -1 } );
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
        dsd_node res = { node_type::xor_, var_index, {} };
        res.fanin.push_back( { 0, -1 } );
        return res;
      }
    }

    return { node_type::none, var_index, {} };
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
      return { node_type::none, -1, {} };
    }

    if ( !eq01 && !eq02 && !eq03 ) // 00 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt11, tt00 );
      }
      dsd_node res = { node_type::or_, new_index, {} };
      res.fanin.push_back( { 0, var_index1 } );
      res.fanin.push_back( { 0, var_index2 } );
      return res;
    }
    else if ( !eq01 && !eq12 && !eq13 ) // 01 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt10 );
      }
      dsd_node res = { node_type::and_, new_index, {} };
      res.fanin.push_back( { 1, var_index1 } );
      res.fanin.push_back( { 0, var_index2 } );
      return res;
    }
    else if ( !eq02 && !eq12 && !eq23 ) // 10 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt10 );
      }
      dsd_node res = { node_type::or_, new_index, {} };
      res.fanin.push_back( { 1, var_index1 } );
      res.fanin.push_back( { 0, var_index2 } );
      return res;
    }
    else if ( !eq03 && !eq13 && !eq23 ) // 11 is different
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt11, tt00 );
      }
      dsd_node res = { node_type::and_, new_index, {} };
      res.fanin.push_back( { 0, var_index1 } );
      res.fanin.push_back( { 0, var_index2 } );
      return res;
    }
    else if ( allow_xor ) // XOR
    {
      if ( func )
      {
        *func = kitty::mux_var( var_index1, tt01, tt00 );
      }
      dsd_node res = { node_type::xor_, new_index, {} };
      res.fanin.push_back( { 0, var_index1 } );
      res.fanin.push_back( { 0, var_index2 } );
      return res;
    }

    return { node_type::none, -1, {} };
  }

  template<class TT>
  int find_unate_var( const TT tt )
  {
    int index;
    for ( index = tt.num_vars() - 1; index > 0; index-- )
    {
      const auto tt0 = kitty::cofactor0( tt, index );
      const auto tt1 = kitty::cofactor1( tt, index );
      if ( ( ( tt0 & tt1 ) == tt0 ) && ( ( tt0 & tt1 ) == tt1 ) )
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

    dsd_node res = { node_type::mux_, index, {} };
    res.fanin.push_back( { 0, index } );

    if ( func0 && func1 )
    {
      *func0 = tt0;
      *func1 = tt1;
    }

    return res;
  }

  /*! \brief Get depth of rule starting from a specific dsd_node.
   *
   *  \param rule rule
   *  \param n dsd_node to start from
   * Returns depth of rule starting from n.
   */
  uint32_t get_depth( rule rule, dsd_node n )
  {
    if ( n.type == node_type::pi_ || n.type == node_type::zero_ )
    {
      return 0;
    }
    uint32_t max_depth;
    uint32_t left_depth = get_depth( rule, rule[n.fanin[0].index] );
    uint32_t right_depth = get_depth( rule, rule[n.fanin[1].index] );
    max_depth = ( left_depth > right_depth ) ? left_depth : right_depth;
    return max_depth + 1;
  }

#pragma region Report
  std::string to_string( node_type t )
  {
    if ( t == node_type::and_ )
      return "*";
    if ( t == node_type::or_ )
      return "+";
    if ( t == node_type::mux_ )
      return "+";
    if ( t == node_type::xor_ )
      return "xor";
    if ( t == node_type::pi_ )
      return "pi";
    if ( t == node_type::none )
      return "none";
    if ( t == node_type::zero_ )
      return "zero";
  }

  void print_dsd_node( dsd_node& n )
  {
    std::cout << n.index << " " << to_string( n.type ) << " ";
    for ( auto elem : n.fanin )
      std::cout << "{" << elem.index << ", " << elem.inv << "}";
    std::cout << "\n";
  }

  void print_rule( rule& r )
  {
    for ( auto elem : r )
      print_dsd_node( elem );
  }

  /*! \brief Print expression of a rule.
   *
   *  \param rule rule.
   *  \param n dsd_node to start from.
   */
  void print_rule( rule rule, dsd_node n )
  {
    if ( n.type == node_type::pi_ )
    {
      std::cout << char( 'a' + n.index );
      return;
    }
    if ( n.type == node_type::zero_ )
    {
      std::cout << "0";
      return;
    }
    else
    {
      std::cout << "(";
      if ( n.fanin[0].inv )
      {
        std::cout << "!";
      }
      if ( n.type == node_type::mux_ )
      {
        std::cout << "!" << char( 'a' + n.fanin[2].index ) << " * ";
      }
      print_rule( rule, rule[n.fanin[0].index] );
      std::cout << " " << to_string( n.type ) << " ";
      if ( n.type == node_type::mux_ )
      {
        std::cout << char( 'a' + n.fanin[2].index ) << " * ";
      }
      if ( n.fanin[1].inv )
      {
        std::cout << "!";
      }
      print_rule( rule, rule[n.fanin[1].index] );
      std::cout << ")";
    }
  }
#pragma endregion

private:
  bool gate_disjoint{ false }; /* flag for gate support*/
  bool not_simmetry{ true };   /* flag for simmetry support*/
  bool permutation{ false };   /* flag for input permutations support */

  std::vector<gate> const& _gates; /* collection of gates */
  composed_list_t _supergates;     /* list of composed_gates */
  lib_rule _dsd_map;               /* hash map for DSD decomposition of gates */
  lib_table _and_table;            /* AND table */
  map_label_gate _label_to_gate;   /* map label to gate */
};

} // namespace mockturtle

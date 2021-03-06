/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "../addons/keras/src/model.h"
#include "./custom.h"
#include "../modules/PhysiCell_settings.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
 
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
	
	// Parameter<double> temp; 
	
	int initial_tumor_radius = parameters.doubles( "tumor_radius" ); 
	
	Cell* pCell; 
	
    double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.8 * 2.0 * cell_radius; 
	
	std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius,initial_tumor_radius);

	std::cout << "Creating cells" << std::endl;
    
    for( int i=0; i < positions.size(); i++ )
    {
        pCell = create_cell(get_cell_definition("default"));
        pCell->functions.volume_update_function = NULL;
        pCell->assign_position( positions[i] );
	}
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 		
	
	return; 
}


std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{  
    return;
} 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


void intracellular_DNN()
{
	static int glc_index = microenvironment.find_density_index( "glucose" );
	static int oxy_index = microenvironment.find_density_index( "oxygen" );
    
	
	#pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        static int biomass_vi = (*all_cells)[i]->custom_data.find_variable_index( "biomass_flux" );
		double cell_volume = (*all_cells)[i]->phenotype.volume.total;
		double glc_val_int = (*all_cells)[i]->nearest_density_vector()[glc_index] * 20;
		double oxy_val_int = (*all_cells)[i]->nearest_density_vector()[oxy_index] * 200;

        float fl_glc = glc_val_int;
        float fl_oxy = oxy_val_int;
        //std::cout << "oxygen --> " << fl_oxy << "        " << "glucose --->    " << fl_glc << std::endl;
		/* auto model = keras2cpp::Model::load("Wild_Type.model"); //model input
		keras2cpp::Tensor in{3}; //
		in.data_ = {fl_glc,0,fl_oxy};
		keras2cpp::Tensor out = model(in); // model evaluation
		//out.print();
		//keras2cpp::Tensor res = out;
        std::vector<double> result; // try to escape copying 
        result = out.result_vector();
        (*all_cells)[i]->custom_data[biomass_vi]  = biomass_creation_flux;
        //std::cout << (*all_cells)[i]->custom_data[biomass_vi] << std::endl; */
        double biomass_creation_flux = 0.035;        
        double volume_increase_ratio = 1 + ( biomass_creation_flux / 60) * 0.01;
        //std::cout << "TOTAL VOLUME OF CELL BEFORE : " << (*all_cells)[i]->phenotype.volume.total << std::endl;
        (*all_cells)[i]->phenotype.volume.multiply_by_ratio(volume_increase_ratio);
        //std::cout << "TOTAL VOLUME OF CELL AFTER : " << (*all_cells)[i]->phenotype.volume.total << std::endl;
        if ( (*all_cells)[i]->phenotype.volume.total > 2494*2)
        {
            //std::cout << "Volume is big enough to divide" << std::endl;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(0,1) = 9e99;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(1,2) = 9e99;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(2,3) = 9e99;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(3,0) = 9e99;
        }
		else
		{
			(*all_cells)[i]->phenotype.cycle.data.transition_rate(0,1) = 0;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(1,2) = 0;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(2,3) = 0;
            (*all_cells)[i]->phenotype.cycle.data.transition_rate(3,0) = 0;
		}
	}
	return;
}
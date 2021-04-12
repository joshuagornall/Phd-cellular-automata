#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <random>
#include <vector>

#define INHI 1	// inhibitor
#define POLY 2	// polymer
#define FILL 5	// filler

const double pi = 3.141592653589793;

typedef std::vector<std::vector<int>> CoordsList;
typedef struct PSD {
	std::vector<double> Psize;
	std::vector<double> Cdf;
}PSD;

std::vector<int> genStructure(int Nx, int Ny, int Nz, std::vector<int> & particle_ID, std::vector<double> & target_PVC, std::vector<PSD> & psd, std::vector<double> & aspect_ratio, double space_step);
void genCylinder(double p_size, double aspect_ratio, int & new_box_size, CoordsList & particle, std::mt19937 & rnd_generator, std::uniform_real_distribution<double> & random_angles, double space_step);
void convertCoordinates(CoordsList & c_particle, CoordsList particle, int p, int Nx, int Ny, int Nz);
void save_coat(std::string filename, std::vector<int> & M, int Nx, int Ny, int Nz);
void read_psd(std::vector<std::string> filename, std::vector<PSD> & psd);
bool check_overlapping(CoordsList & voxel_list, std::vector<int> & structure, int Nx, int Ny, int Nz);
void save_vtk(std::string filename, std::vector<int> & M, int Nx, int Ny, int Nz);

int main(int argc, char **argv) {
	
	//////////////////////////////////////////////////////////////////////////////////
	// SETUP
	
	// geometry
	double space_step = 1.0;	// [microns]
	int Nx = 200;	// [number of voxels]
	int Ny = 200;	// [number of voxels]
	int Nz = 30;	// [number of voxels] - this is the thickness of the coating
	
	// unique ID for each pigment type
	std::vector<int> particle_ID = {INHI};
	
	// volume concentration for each pigment type; follows the order of particle_ID; numbers between 0 and 100
	std::vector<double> particle_PVC = {7.5};
	
	// aspect ratio for each pigment type; follows the order of particle_ID
	std::vector<double> aspect_ratio = {5.0};
	
	// paths to the PSD files (one for each component in the formulation); follows the order of particle_ID
	std::vector<std::string> psd_files = {""};
	
	// name of the output file(s)
	std::string filename_coat = "";
	
	// save coating in VTK format for visualisation
	bool save_vtk_file = true;
	
	//////////////////////////////////////////////////////////////////////////////////
	
	Nx = Nx/space_step;
	Ny = Ny/space_step;
	Nz = Nz/space_step;
	std::vector<PSD> psd(psd_files.size());
	read_psd(psd_files, psd);
	std::vector<int> M(Nx*Ny*Nz);
	M = genStructure(Nx, Ny, Nz, particle_ID, particle_PVC, psd, aspect_ratio, space_step);
	for (int i = 0; i < (Nx*Ny*Nz); i++){
		if (M[i] <= 0) M[i] = POLY;
	}
	save_coat(filename_coat, M, Nx, Ny, Nz);
	if (save_vtk_file){
		save_vtk(filename_coat, M, Nx, Ny, Nz);
	}
	
	return 0;
}

std::vector<int> genStructure(int Nx, int Ny, int Nz, std::vector<int> & particle_ID, std::vector<double> & target_PVC, std::vector<PSD> & psd, std::vector<double> & aspect_ratio, double space_step) {
	
	std::vector<int> structure(Nx*Ny*Nz);
	CoordsList particle;
	CoordsList c_particle;
	
	int rand_loc = 0;
	int temp = 0;
	int num_try = 0;
	int max_num_try = 1000;
	bool overlap_flag = false;
	int new_box_size = Nz + 1;
	double p = 0.0;
	std::vector<double>::iterator first_iter;
	double p_size = 0.0;
	
	std::vector<bool> flag_list((int)(particle_ID.size()), false);
	std::vector<double> current_PVC((int)(particle_ID.size()), 0.0);
	
	std::mt19937 rnd_generator(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<double> random_angles(0.0, pi*2.0000000000001);
	std::uniform_real_distribution<double> random_num(0.0, 1.0000000001);
	
	bool PVC_flag = false;
	
	while (!PVC_flag) {
		for (unsigned k = 0; k < target_PVC.size(); ++k) {
			if (!flag_list[k]) {
				
				p = random_num(rnd_generator);
				
				first_iter = std::find_if(psd[k].Cdf.begin(), psd[k].Cdf.end(), [&p](double c) {return c >= p; });
				
				temp = std::distance(psd[k].Cdf.begin(), first_iter);
				
				p_size = (psd[k].Psize[temp] + psd[k].Psize[temp-1]) / 2.0;
				
				genCylinder(p_size, aspect_ratio[k], new_box_size, particle, rnd_generator, random_angles, space_step);
				
				if (new_box_size < Nz) {
					
					rand_loc = int(random_num(rnd_generator) * (Nx*Ny*(Nz - new_box_size + 1) - 1));
					
					convertCoordinates(c_particle, particle, rand_loc, Nx, Ny, Nz);
					
					overlap_flag = check_overlapping(c_particle, structure, Nx, Ny, Nz);
					
					if (overlap_flag) {
						num_try = 0;
						while (overlap_flag && num_try < max_num_try) {
							rand_loc = int(random_num(rnd_generator) * (Nx*Ny*(Nz - new_box_size + 1) - 1));
							convertCoordinates(c_particle, particle, rand_loc, Nx, Ny, Nz);
							overlap_flag = check_overlapping(c_particle, structure, Nx, Ny, Nz);
							num_try++;
						}
					}
					
					if (!overlap_flag && num_try < max_num_try) {
						for (auto voxel : c_particle) {
							structure[voxel[0] + Nx * voxel[1] + Nx * Ny*voxel[2]] = particle_ID[k];
						}
					}
				
				}
				
				current_PVC[k] = 100.0*(1.0*std::count(structure.begin(), structure.end(), particle_ID[k])) / (1.0*Nx*Ny*Nz);
				
				if ((target_PVC[k] - current_PVC[k]) < 0.0001) flag_list[k] = true;
				
			}
		}
		if (std::all_of(flag_list.begin(), flag_list.end(), [](bool c) {return c; })) PVC_flag = true;
	}
	
	return structure;
}

void genCylinder(double p_size, double aspect_ratio, int & new_box_size, CoordsList & particle, std::mt19937 & rnd_generator, std::uniform_real_distribution<double> & random_angles, double space_step){
	
	double L = p_size/std::sqrt(1.0 + 1.0/(aspect_ratio*aspect_ratio));
	double d = (1.0*L)/(2.0*aspect_ratio);
	if (d < space_step){
		d = space_step;
	}
	if (L < space_step){
		L = space_step;
	}
	L /= space_step;
	d /= space_step;
	
	int box_side = std::round(p_size/space_step);
	
	if (box_side == 0){
		box_side = 1;
	}
	
	double C = box_side/2.0;
	std::vector<int> box(box_side*box_side*box_side);
	
	if (box_side == 1){
		new_box_size = box_side;
		particle.resize(1);
		particle[0] = {{0,0,0}};
		return;
	} 
	
	double nx = cos(random_angles(rnd_generator));
	double ny = cos(random_angles(rnd_generator));
	double nz = cos(random_angles(rnd_generator));
	
	nx = nx/(std::sqrt(std::pow(nx,2.0)+std::pow(ny,2.0)+std::pow(nz,2.0)));
	ny = ny/(std::sqrt(std::pow(nx,2.0)+std::pow(ny,2.0)+std::pow(nz,2.0)));
	nz = nz/(std::sqrt(std::pow(nx,2.0)+std::pow(ny,2.0)+std::pow(nz,2.0)));
	
	double X1 = C + (L/2)*nx;
	double Y1 = C + (L/2)*ny;
	double Z1 = C + (L/2)*nz;
	double X2 = C - (L/2)*nx;
	double Y2 = C - (L/2)*ny;
	double Z2 = C - (L/2)*nz;
	
	double a = 0;
	double b = 0;
	double c = 0;
	double t = L*L*d*d;
	
	for (int z = 0; z < box_side; ++z){
		for (int y = 0; y < box_side; ++y){
			for (int x = 0; x < box_side; ++x){
				a = (y-Y1)*(z-Z2)-(z-Z1)*(y-Y2);
				b = (x-X1)*(z-Z2)-(z-Z1)*(x-X2);
				c = (x-X1)*(y-Y2)-(y-Y1)*(x-X2);
				a = a*a;
				b = b*b;
				c = c*c;
				if (a+b+c<=t && ((x-X1)*nx+(y-Y1)*ny+(z-Z1)*nz)<0 && ((x-X2)*(-nx)+(y-Y2)*(-ny)+(z-Z2)*(-nz))<0){
					box[x + box_side*y + box_side*box_side*z] = 1;
				}
			}
		}
	}
	
	int Zmin = 0;
	while (!std::any_of(box.begin() + box_side * box_side*Zmin, box.begin() + box_side * box_side*(Zmin + 1), [](int c) {return c == 1; })) {
		Zmin++;
	}
	box.erase(box.begin(), box.begin() + (box_side*box_side*(Zmin)));
	int Zmax = box.size() / box_side / box_side;
	while (!std::any_of(box.begin() + box_side * box_side*(Zmax - 1), box.begin() + box_side * box_side*Zmax, [](int c) {return c == 1; })) {
		Zmax--;
	}
	box.erase(box.begin() + (box_side*box_side*(Zmax)), box.end());
	
	new_box_size = box.size() / box_side / box_side;
	
	particle.resize(0);
	int x, y, z;
	for (unsigned k = 0; k < box.size(); ++k) {
		if (box[k] == 1) {
			x = k % box_side;
			y = (k / box_side) % box_side;
			z = k / box_side / box_side;
			particle.push_back({ x, y, z });
		}
	}
	
}

void convertCoordinates(CoordsList & c_particle, CoordsList particle, int p, int Nx, int Ny, int Nz) {
	c_particle = particle;
	for (auto & k : c_particle) {
		k[0] = (k[0] + (p%Nx)) % Nx;
		k[1] = (k[1] + ((p / Nx) % Ny)) % Ny;
		k[2] = (k[2] + (p / Nx / Ny)) % Nz;
	}
}

void save_coat(std::string filename, std::vector<int> & M, int Nx, int Ny, int Nz){
	std::ofstream output(filename + ".coat",std::ios::out | std::ios::binary);
	output.write(reinterpret_cast<char*>(&Nx),sizeof(int));
	output.write(reinterpret_cast<char*>(&Ny),sizeof(int));
	output.write(reinterpret_cast<char*>(&Nz),sizeof(int));
	int state = M[0];
	output.write(reinterpret_cast<char*>(&state),sizeof(int));
	int k = 0;
	output.write(reinterpret_cast<char*>(&k),sizeof(int));
	int end = 0;
	for (k = 1; k < Nx*Ny*Nz; ++k){
		if (M[k] != state){
			end = k-1;
			output.write(reinterpret_cast<char*>(&end),sizeof(int));
			output.write(reinterpret_cast<char*>(&M[k]),sizeof(int));
			output.write(reinterpret_cast<char*>(&k),sizeof(int));
			state = M[k];
		}
	}
	end = Nx*Ny*Nz - 1;
	output.write(reinterpret_cast<char*>(&end),sizeof(int));
}

void read_psd(std::vector<std::string> filename, std::vector<PSD> & psd) {
	double temp_size = 0.0, temp_value = 0.0;
	for (unsigned n = 0; n < filename.size(); ++n) {
		std::ifstream input(filename[n]);
		while (input >> temp_size >> temp_value) {
			psd[n].Psize.push_back(temp_size);
			psd[n].Cdf.push_back(temp_value);
		}
	}
}

bool check_overlapping(CoordsList & voxel_list, std::vector<int> & structure, int Nx, int Ny, int Nz) {
	for (auto voxel : voxel_list) {
		if (structure[voxel[0] + Nx * voxel[1] + Nx * Ny*voxel[2]] != 0) {
			return true;
		}
	}
	return false;
}

void save_vtk(std::string filename, std::vector<int> & M, int Nx, int Ny, int Nz){
	
	std::ofstream output_file(filename + ".vtk",std::ios::out | std::ios::binary);
	
	output_file << "# vtk DataFile Version 2.0" << std::endl;
	output_file << "Coating microstructure" << std::endl;
	output_file << "ASCII" << std::endl;
	output_file << "DATASET STRUCTURED_POINTS" << std::endl;
	output_file << "DIMENSIONS " << Nx+1 << " " << Ny+1 << " " << Nz+1 << std::endl;
	output_file << "ORIGIN 0.0 0.0 0.0" << std::endl;
	output_file << "SPACING 1.0 1.0 1.0" << std::endl;
	output_file << "CELL_DATA " << Nx*Ny*Nz << std::endl;
	output_file << "SCALARS CAstates int 1" << std::endl;
	output_file << "LOOKUP_TABLE default" << std::endl;
	for (int k = 0; k < Nz; ++k){
		for (int j = 0; j < Ny; ++j){
			for (int i = 0; i < Nx; ++i){
				output_file << M[i+Nx*j+Nx*Ny*k] << std::endl;
			}
		}
	}
	
}

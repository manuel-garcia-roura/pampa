#pragma once

#include "NuclearData.hxx"
#include "ConstantNuclearData.hxx"
#include "FeedbackNuclearData.hxx"
#include "PrecursorData.hxx"
#include "ThermalProperties.hxx"
#include "ConstantProperties.hxx"
#include "GraphiteProperties.hxx"
#include "GraphiteMatrixProperties.hxx"
#include "input.hxx"
#include "utils.hxx"

/* The Material class: */
class Material {
   
   private:
      
      /* Nuclear data: */
      NuclearData* nuclear_data = nullptr;
      
      /* Precursor data: */
      PrecursorData* precursor_data = nullptr;
      
      /* Thermal properties: */
      ThermalProperties* thermal_properties = nullptr;
      
      /* Switch for fuel materials: */
      bool fuel = false;
      
      /* Switch for boundary-condition materials: */
      bool bc = false;
      
      /* Switch for split materials: */
      bool split = false;
      
      /* List of submaterials: */
      Array1D<Material*> sub_materials;
      Array1D<int> sub_material_indices;
      
      /* Parent material: */
      Material* parent_material = nullptr;
      int parent_material_index = -1;
   
   public:
      
      /* Material name: */
      const std::string name;
      
      /* The Material constructor: */
      Material(const std::string& name) : name(name) {}
      
      /* The Material copy constructor: */
      Material(const Material& mat, const std::string& name) : nuclear_data(mat.nuclear_data), 
         precursor_data(mat.precursor_data), thermal_properties(mat.thermal_properties), 
         fuel(mat.fuel), bc(mat.bc), split(mat.split), name(name) {}
      
      /* The Material destructor: */
      ~Material() {
         
         /* Free all the material properties: */
         utils::free(&nuclear_data);
         utils::free(&precursor_data);
         utils::free(&thermal_properties);
         
      }
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(const std::string& filename);
      
      /* Read the material from a plain-text input file: */
      int WARN_UNUSED read(std::ifstream& file);
      
      /* Set whether this is a split material: */
      void setSplit(bool split) {this->split = split;}
      
      /* Add a submaterial: */
      void addSubMat(Material* mat, int index) 
         {sub_materials.pushBack(mat); sub_material_indices.pushBack(index);}
      
      /* Set the parent material: */
      void setParentMat(Material* mat, int index) 
         {parent_material = mat; parent_material_index = index;}
      
      /* Check if the material has nuclear data: */
      bool hasNuclearData() const {return nuclear_data != nullptr;}
      
      /* Check if the material has precursor data: */
      bool hasPrecursorData() const {return precursor_data != nullptr;}
      
      /* Check if the material has thermal properties: */
      bool hasThermalProperties() const {return thermal_properties != nullptr;}
      
      /* Check if the material has constant thermal properties: */
      bool hasConstantThermalProperties() const {return thermal_properties->constant;}
      
      /* Check if this is a fuel material: */
      bool isFuel() const {return fuel;}
      
      /* Check if this is a boundary-condition material: */
      bool isBC() const {return bc;}
      
      /* Check if this is a split material: */
      bool isSplit() const {return split;}
      
      /* Get the number of submaterials: */
      int getNumSubMats() const {return sub_materials.size();}
      
      /* Get the submaterials: */
      const Array1D<Material*>& getSubMats() const {return sub_materials;}
      
      /* Get the submaterial indices: */
      const Array1D<int>& getSubMatIndices() const {return sub_material_indices;}
      
      /* Get the parent material: */
      const Material* getParentMat() const {return parent_material;}
      
      /* Get the parent material index: */
      int getParentMatIndex() const {return parent_material_index;}
      
      /* Check the nuclear data to use it in a solver: */
      int WARN_UNUSED checkNuclearData(int num_energy_groups, bool diffusion, bool transient) const 
         {return nuclear_data->check(num_energy_groups, diffusion, transient);}
      
      /* Check the precursor data to use it in a solver: */
      int WARN_UNUSED checkPrecursorData(int num_precursor_groups) const 
         {return precursor_data->check(num_precursor_groups);}
      
      /* Get a total cross section: */
      double sigmaTotal(int g, double T) const {return nuclear_data->sigmaTotal(g, T);}
      
      /* Get a nu-fission cross section: */
      double sigmaNuFission(int g, double T) const {return nuclear_data->sigmaNuFission(g, T);}
      
      /* Get a kappa-fission cross section: */
      double sigmaKappaFission(int g, double T) const 
         {return nuclear_data->sigmaKappaFission(g, T);}
      
      /* Get a transport cross section: */
      double sigmaTransport(int g, double T) const {return nuclear_data->sigmaTransport(g, T);}
      
      /* Get a scattering cross section: */
      double sigmaScattering(int g, int g2, double T) const 
         {return nuclear_data->sigmaScattering(g, g2, T);}
      
      /* Get a diffusion coefficient: */
      double diffusionCoefficient(int g, double T) const 
         {return nuclear_data->diffusionCoefficient(g, T);}
      
      /* Get a prompt-fission spectrum point: */
      double chiPrompt(int g, double T) const {return nuclear_data->chiPrompt(g, T);}
      
      /* Get a delayed-fission spectrum point: */
      double chiDelayed(int g, double T) const {return nuclear_data->chiDelayed(g, T);}
      
      /* Get an effective-fission spectrum point: */
      double chiEffective(int g, double T) const {return nuclear_data->chiEffective(g, T);}
      
      /* Get a neutron velocity: */
      double neutronVelocity(int g, double T) const {return nuclear_data->neutronVelocity(g, T);}
      
      /* Get a precursor decay constant: */
      double lambda(int g) const {return precursor_data->lambda(g);}
      
      /* Get a precursor fraction: */
      double beta(int g) const {return precursor_data->beta(g);}
      
      /* Get the total precursor fraction: */
      double beta() const {return hasPrecursorData() ? precursor_data->beta_total : 0.0;}
      
      /* Get the thermal conductivity: */
      double k(double T) const {return thermal_properties->k(T);}
      
      /* Get the density: */
      double rho(double T) const {return thermal_properties->rho(T);}
      
      /* Get the specific heat capacity: */
      double cp(double T) const {return thermal_properties->cp(T);}
   
};

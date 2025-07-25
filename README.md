# 🌊 Lattice Boltzmann Method Simulation

<div align="center">

![Simulation Banner](https://img.shields.io/badge/Physics-Computational_Fluid_Dynamics-blue?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-3.x-green?style=for-the-badge&logo=python)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)
![Author](https://img.shields.io/badge/Author-Mohamed_Ayoub_Essalami-orange?style=for-the-badge)

**Advanced fluid dynamics simulation demonstrating Kármán vortex street formation behind a cylinder**

</div>

## 📖 Overview

This project implements a comprehensive **Lattice Boltzmann Method (LBM)** simulation that models isothermal fluid flow past a circular cylinder. Using the D2Q9 lattice model, it captures complex fluid phenomena including vortex shedding, boundary layer formation, and turbulent wake development.

The simulation generates a complete visualization suite with **5 different analysis perspectives**, producing both individual and multi-panel videos for detailed flow analysis.

## ✨ Key Features

### 🔬 **Advanced Physics Simulation**
- **D2Q9 Lattice Model**: 2D simulation with 9 velocity directions per lattice point
- **BGK Collision Operator**: Bhatnagar-Gross-Krook approximation for particle collisions
- **Bounce-back Boundaries**: No-slip boundary conditions on cylinder surface
- **Real-time Flow Analysis**: Captures transient and steady-state phenomena

### 📊 **Comprehensive Visualizations**
| Visualization Type | Description | Output Format |
|-------------------|-------------|---------------|
| **Vorticity Field** | Kármán vortex street formation | GIF + MP4 |
| **Velocity Magnitude** | Speed distribution dynamics | GIF + MP4 |
| **Velocity Vectors** | Flow direction field | GIF + MP4 |
| **Density Field** | Pressure variations | GIF + MP4 |
| **Streamlines** | Flow pattern analysis | GIF + MP4 |

### 🎬 **Multi-Format Output**
- **PNG Frames**: Individual time-step images (400+ frames)
- **GIF Animations**: Optimized for web viewing and presentations
- **MP4 Videos**: High-quality videos with title sequences
- **Multi-Panel Video**: Comprehensive analysis combining all visualizations
- **Summary Plot**: Final state analysis with 6 detailed subplots

## 🚀 Quick Start

### Prerequisites
```bash
pip install numpy matplotlib imageio opencv-python pillow tqdm
```

### Running the Simulation
```bash
python "Lattice Boltzmann simulation.py"
```

The simulation will automatically create output directories and generate all visualizations.

## 📁 Output Structure

```
project/
├── frames/              # Individual PNG frames
│   ├── vorticity_*.png
│   ├── velocity_mag_*.png
│   ├── velocity_field_*.png
│   ├── density_*.png
│   └── streamlines_*.png
├── gifs/                # GIF animations
│   ├── 1_vorticity_animation.gif
│   ├── 2_velocity_magnitude_animation.gif
│   ├── 3_velocity_field_animation.gif
│   ├── 4_density_animation.gif
│   └── 5_streamlines_animation.gif
├── videos/              # MP4 videos
│   ├── 1_vorticity_video.mp4
│   ├── 2_velocity_magnitude_video.mp4
│   ├── 3_velocity_field_video.mp4
│   ├── 4_density_video.mp4
│   ├── 5_streamlines_video.mp4
│   └── comprehensive_multi_panel.mp4
└── LBM_Comprehensive_Summary.png
```

## 🎯 Visualization Gallery

### 🌀 **Vorticity Field - Kármán Vortex Street**
*Captures the alternating vortex shedding pattern characteristic of flow past bluff bodies*

**Key Insights:**
- Periodic vortex formation and shedding
- Alternating positive/negative vorticity regions
- Classic Kármán vortex street pattern

### 🏃‍♂️ **Velocity Magnitude Distribution**
*Shows speed variations across the flow field*

**Key Insights:**
- High-speed regions around cylinder edges
- Wake region with reduced velocities
- Acceleration zones and flow separation points

### 🧭 **Velocity Vector Field**
*Displays flow direction and magnitude using arrow vectors*

**Key Insights:**
- Flow streamline deflection around cylinder
- Recirculation zones in the wake
- Boundary layer development

### 🌡️ **Density Field (Pressure Variations)**
*Reveals pressure distribution and fluctuations*

**Key Insights:**
- High-pressure stagnation point
- Low-pressure wake region
- Pressure oscillations due to vortex shedding

### 🌊 **Streamlines Flow Pattern**
*Traces fluid particle paths through the flow field*

**Key Insights:**
- Smooth flow around cylinder front
- Flow separation points
- Complex wake structure

## 🔬 Physics Phenomena Captured

| Phenomenon | Description | Observable in |
|------------|-------------|---------------|
| **Kármán Vortex Street** | Alternating vortex shedding pattern | Vorticity, Streamlines |
| **Boundary Layer Formation** | Viscous effects near surfaces | Velocity Field, Streamlines |
| **Flow Separation** | Detachment from cylinder surface | All visualizations |
| **Wake Turbulence** | Chaotic downstream flow | Vorticity, Velocity |
| **Pressure Variations** | Density fluctuations | Density Field |
| **Velocity Gradients** | Speed variations across field | Velocity Magnitude |

## ⚙️ Technical Specifications

### **Simulation Parameters**
- **Grid Resolution**: 400 × 100 lattice points
- **Time Steps**: 4,000 iterations
- **Relaxation Time (τ)**: 0.6
- **Average Density (ρ₀)**: 0.6
- **Reynolds Number**: ~100-200 (depending on flow conditions)

### **Computational Method**
- **Lattice**: D2Q9 (9 discrete velocities in 2D)
- **Collision**: BGK single-relaxation-time model
- **Boundary**: Bounce-back for cylinder surface
- **Streaming**: Standard LBM streaming step
- **Initial Conditions**: Uniform flow with small perturbations

### **Performance Metrics**
- **Lattice Points**: 40,000 total
- **Memory Usage**: ~3 MB for distribution functions
- **Output Files**: 16 total (5 GIFs + 6 videos + 5 frame sets)
- **Processing Time**: ~2-5 minutes (depending on hardware)

## 📚 Recommended Viewing Sequence

1. **🎬 Start with GIFs** → Quick overview of each phenomenon
2. **📹 Individual Videos** → Detailed analysis of specific aspects
3. **🎞️ Multi-Panel Video** → Comprehensive simultaneous comparison
4. **📊 Summary Plot** → Final state quantitative analysis

## 🛠️ Advanced Usage

### **Customizing Parameters**
Modify simulation parameters in the `LBMSimulation.__init__()` method:
```python
self.Nx = 400      # Grid resolution x-direction
self.Ny = 100      # Grid resolution y-direction
self.tau = 0.6     # Relaxation time (affects viscosity)
self.Nt = 4000     # Number of time steps
```

### **Output Customization**
Adjust visualization settings:
- Frame skip rate: Modify `if ((it % 25) == 0)` in `run_simulation()`
- Color maps: Change `cmap` parameters in `save_frame()`
- Video quality: Adjust `dpi` and `fps` parameters

## 🔧 Troubleshooting

### **Common Issues**

| Issue | Solution |
|-------|----------|
| **Import Errors** | Install missing packages: `pip install package_name` |
| **Memory Issues** | Reduce grid size (`Nx`, `Ny`) or time steps (`Nt`) |
| **Slow Performance** | Decrease frame save frequency or reduce resolution |
| **File Access Errors** | Check write permissions and available disk space |

### **System Requirements**
- **Python**: 3.7+
- **RAM**: 4GB minimum (8GB recommended)
- **Storage**: 500MB free space for output files
- **OS**: Windows, macOS, or Linux

## 🤝 Contributing

Contributions are welcome! Areas for improvement:
- Additional boundary conditions
- 3D visualization capabilities
- Interactive parameter adjustment
- Performance optimizations
- Additional physics models

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👨‍💻 Author

**Mohamed Ayoub Essalami**
- Computational Fluid Dynamics Researcher
- Lattice Boltzmann Method Specialist

## 📚 References

- Chen, S., & Doolen, G. D. (1998). Lattice Boltzmann method for fluid flows. *Annual Review of Fluid Mechanics*, 30(1), 329-364.
- Succi, S. (2001). *The lattice Boltzmann equation: for fluid dynamics and beyond*. Oxford University Press.
- Krüger, T., et al. (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.

## 🎯 Acknowledgments

- Original inspiration from Philip Mocz's Princeton University LBM tutorial
- Extended and enhanced with comprehensive visualization suite
- Optimized for educational and research applications

---

<div align="center">

**⭐ If this project helped you, please give it a star! ⭐**

*Built with passion for computational fluid dynamics and scientific visualization*

</div>
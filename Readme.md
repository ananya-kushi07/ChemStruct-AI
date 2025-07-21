# 🚀 ChemStruct AI: Dynamic Molecular Visualization & Reaction Simulation

> **Unlock the Dynamics of Molecules**  
ChemStruct AI is a cutting-edge web application built to help students and researchers intuitively explore molecular structures and their dynamic reactivity. By integrating 3D rendering and Graph Neural Networks (GNN), ChemStruct AI transforms static chemical data into interactive animated simulations.

---

## ✨ Problem Statement

Traditional chemical diagrams and tools often lack the dynamic, interactive feedback necessary for intuitive understanding. ChemStruct AI addresses this gap by offering a real-time visualization and reactivity prediction platform accessible through a web browser.

---

## 💡 Value Proposition

- **Enhanced Learning:** Makes chemistry engaging and accessible with real-time 3D animations.
- **Accelerated Research:** Fast, AI-powered predictions of reactive bonds for researchers.
- **Web Accessibility:** No need for heavy installations—just open in your browser.
- **High Engagement:** Visually stunning animations deliver the "wow" factor in education or demos.

---

## 🌟 Core Features (Hackathon MVP)

- 🔤 Input molecular formulas or **SMILES strings**
- 🔍 **Real-time 3D molecular visualization** with WebGL & Three.js
- 🧠 **AI-Powered Bond Reactivity Prediction** using Graph Neural Networks
- 🎬 **Animated Reaction Simulation** based on prediction
- 🎯 Predicted bond clearly **highlighted** in the 3D model and **displayed textually**
- 📱 **Mobile-responsive** design with a clean modern UI/UX

---

## 🛠️ Technical Stack

### Frontend
- **React.js** – Dynamic UI
- **Three.js / WebGL** – High-performance 3D rendering
- **HTML5 / CSS3** – Modern styling

### Backend
- **Flask (Python)** – REST API & ML orchestration
- **RDKit** – Molecular structure parsing & feature extraction

### Machine Learning
- **PyTorch / PyTorch Geometric** – GNN implementation
- **GNNs** – Reactivity prediction using molecular graphs

### APIs
- **PubChem** *(planned)* – Molecular data enrichment
- **Internal ML API** – Flask-based ML model inference

---

## 🧱 Architecture

```mermaid
graph TD
    A[User] --> B(Frontend - React.js)
    B --> C(Backend - Flask API)
    C --> D(ML Model - PyTorch GNN)
    D --> C
    C --> B
    B --> E[3D Molecular Visualization & Animation]

// src/App.js
import React, { useState } from "react";
import MoleculeInput from "./components/MoleculeInput";
import MolecularViewer from "./components/MolecularViewer";
import { predictBondReactivity } from './api/Api';
import "./App.css";

function App() {
  const [moleculeData, setMoleculeData] = useState(null);
  const [prediction, setPrediction] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  const handleSimulate = async (smiles) => {
    setLoading(true);
    setError("");
    setPrediction(null);
    try {
      const response = await predictBondReactivity(smiles);
      setMoleculeData(response.molecule_info);
      setPrediction(response.prediction);
    } catch (err) {
      console.error(err);
      setError("Prediction failed. Please check the input.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="app-container">
      <header>
        <h1>ChemStruct AI</h1>
        <p>Predict & Visualize Molecular Bond Reactivity</p>
      </header>

      <main className="main-layout">
        <div className="left-panel">
          <MoleculeInput onSubmit={handleSimulate} loading={loading} />
          {error && <div className="error">{error}</div>}

          {moleculeData && prediction && (
            <div className="prediction-panel">
              <h3>Prediction Output</h3>
              <p><strong>Formula:</strong> {moleculeData.formula}</p>
              <p><strong>Weight:</strong> {moleculeData.molecular_weight}</p>
              <p><strong>Atoms:</strong> {moleculeData.num_atoms}</p>
              <p><strong>Bonds:</strong> {moleculeData.num_bonds}</p>
              <p className="highlight">
                🔥 Predicted Break: Atoms {prediction.predicted_bond[0]} - {prediction.predicted_bond[1]}
              </p>
              <p>Confidence: {(prediction.confidence * 100).toFixed(2)}%</p>
            </div>
          )}
        </div>

        <div className="viewer-panel">
          <MolecularViewer prediction={prediction} />
        </div>
      </main>

      <footer>
        <p>💡 Built with React + Three.js + Flask + RDKit</p>
      </footer>
    </div>
  );
}

export default App;

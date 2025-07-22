// src/components/MoleculeInput.js
import React, { useState } from 'react';

function MoleculeInput({ onSubmit, loading }) {
  const [input, setInput] = useState('');
  const sampleMolecules = ['ethanol', 'CCO', 'C2H6O', 'H2O'];

  const handleClick = () => {
    if (!input.trim()) return;
    if (onSubmit) onSubmit(input); // send only input to App.js
  };

  return (
    <div style={{ marginBottom: '1rem' }}>
      <input
        type="text"
        placeholder="Enter SMILES or molecule name..."
        value={input}
        onChange={(e) => setInput(e.target.value)}
        disabled={loading}
        style={{ padding: '8px', width: '60%' }}
      />
      <button onClick={handleClick} disabled={loading} style={{ marginLeft: '10px' }}>
        {loading ? 'Predicting...' : 'Simulate'}
      </button>

      <div style={{ marginTop: '10px' }}>
        <strong>Sample molecules:</strong>
        {sampleMolecules.map((mol, idx) => (
          <button
            key={idx}
            onClick={() => setInput(mol)}
            disabled={loading}
            style={{ marginLeft: '8px' }}
          >
            {mol}
          </button>
        ))}
      </div>
    </div>
  );
}

export default MoleculeInput;

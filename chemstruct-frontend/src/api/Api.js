// src/api/Api.js

// ✅ Make sure URL is lowercase '/api'
const BASE_URL = "http://127.0.0.1:5000/api";

export async function predictBondReactivity(smiles) {
  try {
    const response = await fetch(`${BASE_URL}/predict_bond_reactivity`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ smiles: smiles }), // frontend sends JSON
    });

    if (!response.ok) {
      const err = await response.text();
      throw new Error(`Backend error: ${err}`);
    }

    // parse JSON
    const data = await response.json();
    return data;

  } catch (error) {
    console.error("❌ Error while calling backend:", error);
    throw error;
  }
}

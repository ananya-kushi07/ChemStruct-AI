// src/components/MolecularViewer.js
import React, { useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';

// 🔹 Bond component
function Bond({ from, to, isHighlighted }) {
  const meshRef = useRef();
  const start = new THREE.Vector3(...from);
  const end = new THREE.Vector3(...to);
  const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  const dir = new THREE.Vector3().subVectors(end, start);
  const len = dir.length();

  const quaternion = new THREE.Quaternion();
  quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());

  useFrame((state) => {
    if (isHighlighted && meshRef.current) {
      const time = state.clock.getElapsedTime();
      const pulse = 0.08 + 0.02 * Math.sin(time * 5);
      meshRef.current.scale.set(pulse / 0.08, 1, pulse / 0.08);
    }
  });

  return (
    <mesh ref={meshRef} position={mid} quaternion={quaternion}>
      <cylinderGeometry args={[0.08, 0.08, len, 16]} />
      <meshStandardMaterial
        color={isHighlighted ? 'orange' : 'gray'}
        emissive={isHighlighted ? 'yellow' : 'black'}
      />
    </mesh>
  );
}

// 🔹 Main viewer
function MolecularViewer({ prediction }) {
  // For now, use dummy test bonds
  const bonds = [
    { from: [0, 0, 0], to: [0, 1, 0] },
    { from: [0, 1, 0], to: [1, 1, 0] },
  ];

  // Decide highlight based on prediction
  const highlightIndices = prediction?.predicted_bond || [];

  return (
    <div style={{ width: '100%', height: '400px' }}>
      <Canvas camera={{ position: [2, 2, 4] }}>
        <ambientLight intensity={0.5} />
        <directionalLight position={[5, 5, 5]} intensity={1} />

        {bonds.map((b, idx) => (
          <Bond
            key={idx}
            from={b.from}
            to={b.to}
            isHighlighted={highlightIndices.includes(idx)}
          />
        ))}

        <OrbitControls enableZoom enableRotate />
      </Canvas>
    </div>
  );
}

export default MolecularViewer;

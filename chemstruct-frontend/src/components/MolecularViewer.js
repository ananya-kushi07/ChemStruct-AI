// src/components/MolecularViewer.js
import React, { useRef, useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import * as THREE from 'three';

// atom color map
const atomColors = {
  C: '#4F4F4F',  // dark gray
  H: '#FFFFFF',  // white
  O: '#FF0000',  // red
  N: '#0000FF',  // blue
  default: '#AAAAAA'
};

// 🌟 Atom component
function Atom({ position, element }) {
  const color = atomColors[element] || atomColors.default;
  return (
    <mesh position={position}>
      <sphereGeometry args={[0.2, 32, 32]} />
      <meshStandardMaterial color={color} roughness={0.2} metalness={0.5} />
    </mesh>
  );
}

// 🌟 Bond component
function Bond({ from, to, isHighlighted }) {
  const meshRef = useRef();
  const start = new THREE.Vector3(...from);
  const end = new THREE.Vector3(...to);
  const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);
  const dir = new THREE.Vector3().subVectors(end, start);
  const len = dir.length();
  const quaternion = new THREE.Quaternion();
  quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());

  // pulse animation
  useFrame((state) => {
    if (isHighlighted && meshRef.current) {
      const t = state.clock.getElapsedTime();
      const pulse = 0.08 + 0.02 * Math.sin(t * 5);
      meshRef.current.scale.set(pulse / 0.08, 1, pulse / 0.08);
    }
  });

  return (
    <mesh ref={meshRef} position={mid} quaternion={quaternion}>
      <cylinderGeometry args={[0.08, 0.08, len, 16]} />
      <meshStandardMaterial
        color={isHighlighted ? 'orange' : '#888'}
        emissive={isHighlighted ? 'yellow' : 'black'}
        roughness={0.3}
        metalness={0.7}
      />
    </mesh>
  );
}

// 🌟 Main Viewer
function MolecularViewer({ prediction }) {
  // Dummy atoms/bonds – replace with real from backend later
  const atoms = [
    { pos: [0, 0, 0], el: 'C' },
    { pos: [1, 0, 0], el: 'O' },
    { pos: [-1, 0, 0], el: 'H' },
  ];
  const bonds = [
    { from: [0, 0, 0], to: [1, 0, 0], isHighlighted: prediction ? prediction.predicted_bond.includes(0) : false },
    { from: [0, 0, 0], to: [-1, 0, 0], isHighlighted: false },
  ];

  // UI states
  const [autoRotate, setAutoRotate] = useState(true);
  const [showBonds, setShowBonds] = useState(true);

  return (
    <div style={{ position: 'relative', width: '100%', height: '400px' }}>
      {/* Fancy UI controls */}
      <div style={{
        position: 'absolute', top: 10, right: 10, zIndex: 10,
        background: 'rgba(255,255,255,0.8)', padding: '6px 10px', borderRadius: '8px'
      }}>
        <button onClick={() => setAutoRotate(!autoRotate)} style={{ marginRight: '6px' }}>
          {autoRotate ? '⏸ Stop Rotate' : '🔄 Auto Rotate'}
        </button>
        <button onClick={() => setShowBonds(!showBonds)}>
          {showBonds ? 'Hide Bonds' : 'Show Bonds'}
        </button>
      </div>

      <Canvas camera={{ position: [2, 2, 4], fov: 50 }}>
        <ambientLight intensity={0.4} />
        <directionalLight position={[5, 5, 5]} intensity={1} castShadow />
        <Environment preset="city" /> {/* ✅ environment preset */}

        {atoms.map((a, i) => <Atom key={i} position={a.pos} element={a.el} />)}
        {showBonds && bonds.map((b, i) => (
          <Bond key={i} from={b.from} to={b.to} isHighlighted={b.isHighlighted} />
        ))}

        <OrbitControls enableZoom enablePan autoRotate={autoRotate} />
      </Canvas>
    </div>
  );
}

export default MolecularViewer;

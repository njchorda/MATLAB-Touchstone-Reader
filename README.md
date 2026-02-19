# SPARAMS — MATLAB S-Parameter Toolbox

A MATLAB class for reading, manipulating, converting, and plotting S-parameter data from Touchstone files (`.s1p`, `.s2p`, `.s3p`, `.s4p`). Supports up to 4-port measurements in **RI**, **MA**, and **dB/Angle** formats.

---

## Getting Started

### Load from a Touchstone file
```matlab
s = SPARAMS('myfile.s2p');
```

### Create an empty object and assign data manually
```matlab
s = SPARAMS();
s.setNumPorts(2);
s.setZ0(50);
s.f   = freq;   % frequency vector in Hz
s.S11 = s11;    % complex S-parameter vectors
s.S21 = s21;
% etc.
```

### Copy an object (no pointer referencing)
```matlab
s2 = s.copyobj();
```

---

## Properties

| Property | Description |
|---|---|
| `f` | Frequency vector (Hz) |
| `S11`, `S12`, ... `S44` | Complex S-parameters (up to 4-port) |
| `Z0` | System reference impedance (e.g. 50 Ω) |
| `numPorts` | Number of ports (1–4) |

---

## Methods

### Setup & Configuration

```matlab
s.setFile('filename.s2p')      % Load a new file into an existing object
s.setZ0(50)                    % Set system impedance manually
s.setNumPorts(2)               % Set number of ports manually
s.setFreqUnits('GHz')          % Set frequency axis units for plotting
                               % Options: 'Hz', 'kHz', 'MHz', 'GHz'
```

### Accessing Data

```matlab
s11 = s.S11;                   % Get S11 as complex (R + jI)
[S, fActual] = s.toMat(freq)   % Get full S-matrix at a specified frequency
                               % (uses nearest frequency if exact match not found)
```

### Plotting

```matlab
s.plotdB('S11')                        % Plot a single parameter in dB
s.plotdB({'S11', 'S21'})               % Plot multiple parameters in dB
s.plotAlldB()                          % Plot all S-parameters in dB
s.plotPolar({'S11', 'S21'})            % Polar plot
s.plotSmith('S11')                     % Smith chart plot
s.plotDispersion()                     % Dispersion diagram (for unit cells)
s.plotDispersion('flip')               % Dispersion with flipped axes
```

### Parameter Conversions

```matlab
[Zin1, Zin2] = s.toInputImpedance()           % Input impedance at each port
[A, B, C, D] = s.toABCDparams()               % ABCD parameters (all frequencies)
[A, B, C, D] = s.toABCDparams(freq)           % ABCD parameters at one frequency
[T11,T12,T21,T22] = s.toTParams()             % T-parameters (all frequencies)
[T11,T12,T21,T22] = s.toTParams(freq)         % T-parameters at one frequency
[Z11,Z12,Z21,Z22] = s.toZparams()             % Z-parameters
[Y11,Y12,Y21,Y22] = s.toYparams()             % Y-parameters
```

### Manipulation

```matlab
s.renorm(Z0new)      % Renormalize to a new reference impedance
                     % (also updates s.Z0 automatically)
s.cascade(N)         % Cascade the unit cell N times (2-port only)
                     % Note: overwrites original S-parameters
                     % Use copyobj() first to preserve originals
```

### Export

```matlab
s.writeCSV('output.csv')    % Export S-parameter data to CSV
s.writeSNP('output')        % Write a new Touchstone file
                            % Extension (.s1p/.s2p/etc.) is added automatically
                            % Format is GHz, RI
```

---

## Static Methods

### De-embedding
Remove error networks from a measurement to isolate the DUT:

```
P1 <-- err1 <--> DUT <--> err2 --> P2
```

```matlab
sDUT = SPARAMS.deembed(err1, measured, err2);
```

Returns a new `SPARAMS` object representing the DUT.

### General Matrix Conversions

```matlab
Z = SPARAMS.s2z(S, Z0)              % S-matrix → Z-matrix
S = SPARAMS.z2s(Z, Z0)              % Z-matrix → S-matrix
[A,B,C,D] = SPARAMS.s2abcd(S, Z0)  % S-matrix → ABCD
[S11,S21,S12,S22] = SPARAMS.abcd2s(A, B, C, D, Z0)  % ABCD → S-parameters
```

### Comparing Two Objects

```matlab
SPARAMS.plotAnydB(obj1, {'S11','S21'}, obj2, {'S11','S21'});
```

---

## Example Workflow

```matlab
% Load a 2-port measurement
s = SPARAMS('dut.s2p');
s.setFreqUnits('GHz');

% Plot S11 and S21
figure;
s.plotdB({'S11', 'S21'});
title('DUT S-Parameters');

% Get ABCD parameters and cascade 3x
sc = s.copyobj();
sc.cascade(3);

% Export result
sc.writeSNP('cascaded_3x');
```

---

## Supported File Formats

| Format | Description |
|---|---|
| `.s1p` | 1-port Touchstone |
| `.s2p` | 2-port Touchstone |
| `.s3p` | 3-port Touchstone |
| `.s4p` | 4-port Touchstone |

Data formats supported within files: **RI** (Real/Imaginary), **MA** (Magnitude/Angle), **dB** (dB/Angle).

---

## Notes

- `cascade()` and `plotDispersion()` are intended for **2-port unit cell** measurements.
- `renorm()` is fully implemented for **1- and 2-port** measurements.
- `toZparams()`, `toYparams()`, `toABCDparams()`, `toInputImpedance()` require **2-port** data.
- When using `cascade()`, it is recommended to `copyobj()` first, as the method overwrites the object's S-parameters.

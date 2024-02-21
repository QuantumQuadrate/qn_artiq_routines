import numpy as np
from scipy.optimize import minimize


# Function to generate random measurement data
def generate_data(num_measurements, num_angles, lp_angle = 4*np.pi/3):
    # Generate random angles for the quarter wave plate and half wave plate
    quarter_wave_angles = np.random.uniform(0, np.pi, num_angles)
    half_wave_angles = np.random.uniform(0, np.pi, num_angles)

    # Generate random measurements for each angle
    measurements = []
    input_states = []
    input_state = np.random.rand(2) + 1j * np.random.rand(2)
    for q_angle, h_angle in zip(quarter_wave_angles, half_wave_angles):
        # Define the Jones matrices for the quarter wave plate and half wave plate
        quarter_wave_matrix = np.exp(complex(-1j*np.pi/4))*np.array([[np.cos(2 * q_angle), np.sin(2 * q_angle)],
                                        [np.sin(2 * q_angle), -np.cos(2 * q_angle)]])
        half_wave_matrix = np.exp(complex(-1j*np.pi/4))*np.array([[np.cos(2 * h_angle), np.sin(2 * h_angle)],
                                     [np.sin(2 * h_angle), -np.cos(2 * h_angle)]])
        linear_polarizer_matrix = np.array([[np.cos(2 * lp_angle), np.sin(2 * lp_angle)],
                                            [np.sin(2 * lp_angle), -np.cos(2 * lp_angle)]])
        # Generate random input states
        input_state /= np.linalg.norm(input_state)  # Normalize to ensure it's a valid quantum state
        input_states.append(input_state)

        # Calculate the predicted output intensities
        predicted_intensity = np.abs(np.dot(half_wave_matrix, np.dot(quarter_wave_matrix, input_state)))[0] ** 2

        # Add noise to the measurements
        measurements.append(predicted_intensity + np.random.normal(scale=0.01))

    return quarter_wave_angles, half_wave_angles, input_states, measurements


# Objective function to minimize
def objective_function(angles, input_states, measurements, linear_polarizer_angle = 4*np.pi/3):
    quarter_wave_angles, half_wave_angles = np.split(angles, 2)
    error = 0

    for q_angle, h_angle, input_state, measured_intensity in zip(quarter_wave_angles, half_wave_angles, input_states,
                                                                 measurements):
        # Define the Jones matrices for the quarter wave plate and half wave plate
        quarter_wave_matrix = np.array([[np.cos(2 * q_angle), np.sin(2 * q_angle)],
                                        [np.sin(2 * q_angle), -np.cos(2 * q_angle)]])
        half_wave_matrix = np.array([[np.cos(2 * h_angle), np.sin(2 * h_angle)],
                                     [np.sin(2 * h_angle), -np.cos(2 * h_angle)]])
        lp_matrix = np.array([[np.cos(2 * linear_polarizer_angle), np.sin(2 * linear_polarizer_angle)],
                                            [np.sin(2 * linear_polarizer_angle), -np.cos(2 * linear_polarizer_angle)]])

        # Calculate the predicted output intensity
        predicted_intensity = np.abs((np.dot(half_wave_matrix, np.dot(quarter_wave_matrix, input_state)))) ** 2

        # Calculate the squared difference between predicted and measured intensity
        error += np.sum((predicted_intensity - measured_intensity) ** 2)

    return error


# Generate random measurement data
num_measurements = 1
num_angles = 5
quarter_wave_angles, half_wave_angles, input_states, measurements = generate_data(num_measurements, num_angles)

# Initial guess for the quarter wave plate and half wave plate angles
initial_guess = np.random.uniform(0, np.pi, 2 * num_angles)

# Minimize the objective function to find the optimal angles
result = minimize(objective_function, initial_guess, args=(input_states, measurements,), method='L-BFGS-B')

# Extract the optimal angles
optimal_angles = result.x
optimal_quarter_wave_angles, optimal_half_wave_angles = np.split(optimal_angles, 2)
def radians_to_degrees(angle_radians):
    angle_degrees = np.degrees(angle_radians)
    return angle_degrees
# Print the optimal quarter wave plate and half wave plate angles
print("Optimal quarter wave plate angles:")
for i, angle in enumerate(optimal_quarter_wave_angles):
    print(f"Angle {i + 1}: {radians_to_degrees(angle)}")

print("Optimal half wave plate angles:")
for i, angle in enumerate(optimal_half_wave_angles):
    print(f"Angle {i + 1}: {radians_to_degrees(angle)}")

def calculate_predicted_output_intensity(quarter_wave_angle_degrees, half_wave_angle_degrees, linear_polarizer_angle_degrees, input_state):
    # Convert angles from degrees to radians
    quarter_wave_angle = np.radians(quarter_wave_angle_degrees)
    half_wave_angle = np.radians(half_wave_angle_degrees)
    linear_polarizer_angle = np.radians(linear_polarizer_angle_degrees)

    # Define Jones matrices for quarter wave plate, half wave plate, and linear polarizer
    quarter_wave_matrix = np.array([[np.cos(2 * quarter_wave_angle), np.sin(2 * quarter_wave_angle)],
                                    [np.sin(2 * quarter_wave_angle), -np.cos(2 * quarter_wave_angle)]])
    half_wave_matrix = np.array([[np.cos(2 * half_wave_angle), np.sin(2 * half_wave_angle)],
                                 [np.sin(2 * half_wave_angle), -np.cos(2 * half_wave_angle)]])
    linear_polarizer_matrix = np.array([[np.cos(2 * linear_polarizer_angle), np.sin(2 * linear_polarizer_angle)],
                                        [np.sin(2 * linear_polarizer_angle), -np.cos(2 * linear_polarizer_angle)]])

    # Calculate the predicted output intensity
    predicted_output_intensity = np.abs(np.dot(linear_polarizer_matrix, np.dot(half_wave_matrix, np.dot(quarter_wave_matrix, input_state))))[0] ** 2

    return predicted_output_intensity

qw_degrees = radians_to_degrees(optimal_quarter_wave_angles)
hw_degrees = radians_to_degrees(optimal_half_wave_angles)

# Define the Jones matrices for the quarter wave plate and half wave plate
predicted_intensity = []
for i, q, h in zip(input_states, optimal_quarter_wave_angles, optimal_half_wave_angles):
    print(q,h)
    quarter_wave_matrix = np.exp(complex(-1j*np.pi/4))*np.array([[np.cos(2 * q), np.sin(2 * q)],
                                        [np.sin(2 * q), -np.cos(2 * q)]])
    half_wave_matrix = np.exp(complex(-1j*np.pi/2))*np.array([[np.cos(2 * h), np.sin(2 * h)],
                            [np.sin(2 * h), -np.cos(2 * h)]])
    predicted_intensity.append(np.abs(np.dot(half_wave_matrix, np.dot(quarter_wave_matrix, i)))[0]**2)

print(predicted_intensity)
print(input_states)

# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt

# load the dataset
dataset = loadtxt('WT_in_silico_data.csv', delimiter=',')
sorted_data = np.array(dataset)
data = np.array(dataset)
np.random.shuffle(data)

percent_training_data = 0.8
splitter = round(data.shape[0] * percent_training_data)

training_data = data[:splitter]
test_data = data[splitter:]


# glu_vals = dataset[:,0]
# glt_vals = dataset[:,2]
# biomass_vals = dataset[:,3]


# X, Y = np.meshgrid(glu_vals, glt_vals)
# biomass_vals3d = np.array([biomass_vals,biomass_vals])

# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(projection='3d')
# ax.scatter(glu_vals, glt_vals, biomass_vals)
# ax.set_xlabel('oxygen')
# ax.set_ylabel('glucose')
# ax.set_zlabel('biomass')
# ax.set_zlim3d(0,7)
# plt.show()



# split into input (X) and output (y) variables
X = training_data[:,0:3]
y = training_data[:,3]


# define the keras model
model = Sequential()
model.add(Dense(4, input_shape=(3,), activation='relu'))
model.add(Dense(12, activation='relu'))
model.add(Dense(1, activation='relu'))
# compile the keras model
model.compile(loss='mse', optimizer='adam', metrics=['mae'])
# fit the keras model on the dataset
history = model.fit(X, y, epochs=150, batch_size=100)

plt.plot(history.history['mae'],'o',color='black')
plt.title('mean absolutue error')
plt.ylabel('mea')
plt.xlabel('number of epochs')
plt.show()



#%%
x_test_data = test_data[:,0:3]
y_test_data = test_data[:,3]


# evaluate the keras model
print('EVALUATION')
_, accuracy = model.evaluate(x_test_data, y_test_data)


#%%
oxygen_values_500 = []
biomass_reals_500 = []
biomass_predictions_500 = []
for i in range(0,20):
    oxygen_value = sorted_data[i,2]
    biomass_real = sorted_data[i,3]
    biomass_predicted = sol = model.predict([[500., 0 ,oxygen_value]])
    print(biomass_predicted)
    oxygen_values_500.append(oxygen_value)
    biomass_reals_500.append(biomass_real)
    biomass_predictions_500.append(biomass_predicted[0])

plt.plot(oxygen_values_500,biomass_reals_500,'ko')
plt.plot(oxygen_values_500,biomass_predictions_500,'b')
plt.xlabel('oxygen lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 500')



oxygen_values_236 = []
biomass_reals_236 = []
biomass_predictions_236 = []
for i in range(4300,4320):
    oxygen_value = sorted_data[i,2]
    biomass_real = sorted_data[i,3]
    biomass_predicted = sol = model.predict([[236.84, 0 ,oxygen_value]])
    print(biomass_predicted)
    oxygen_values_236.append(oxygen_value)
    biomass_reals_236.append(biomass_real)
    biomass_predictions_236.append(biomass_predicted[0])

plt.plot(oxygen_values_236,biomass_reals_236,'ko')
plt.plot(oxygen_values_236,biomass_predictions_236,'b')
plt.xlabel('oxygen lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 236')




oxygen_values_26 = []
biomass_reals_26 = []
biomass_predictions_26 = []
for i in range(7300,7320):
    oxygen_value = sorted_data[i,2]
    biomass_real = sorted_data[i,3]
    biomass_predicted = sol = model.predict([[26.315, 0 ,oxygen_value]])
    print(biomass_predicted)
    oxygen_values_26.append(oxygen_value)
    biomass_reals_26.append(biomass_real)
    biomass_predictions_26.append(biomass_predicted[0])

plt.plot(oxygen_values_26,biomass_reals_26,'ko')
plt.plot(oxygen_values_26,biomass_predictions_26,'b')
plt.xlabel('oxygen lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 26')


oxygen_values_0 = []
biomass_reals_0 = []
biomass_predictions_0 = []
for i in range(7960,7980):
    oxygen_value = sorted_data[i,2]
    biomass_real = sorted_data[i,3]
    biomass_predicted = sol = model.predict([[0, 0 ,oxygen_value]])
    print(biomass_predicted)
    oxygen_values_0.append(oxygen_value)
    biomass_reals_0.append(biomass_real)
    biomass_predictions_0.append(biomass_predicted[0])

plt.plot(oxygen_values_0,biomass_reals_0,'ko')
plt.plot(oxygen_values_0,biomass_predictions_0,'b')
plt.xlabel('oxygen lower boundaries')
plt.ylabel('biomass growth')
plt.title('Glucose lb = 0')





#%%
#save model
from keras2cpp import export_model
export_model(model, 'Wild_Type.model')
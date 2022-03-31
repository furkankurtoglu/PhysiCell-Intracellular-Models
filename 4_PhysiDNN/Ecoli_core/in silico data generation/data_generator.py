# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt

# load the dataset
dataset = loadtxt('in_silico_data.csv', delimiter=',')
data = np.array(dataset)
np.random.shuffle(data)

percent_training_data = 0.8
splitter = round(data.shape[0] * percent_training_data)

training_data = data[:splitter]
test_data = data[splitter:]


glu_vals = dataset[:,0]
glt_vals = dataset[:,1]
biomass_vals = dataset[:,3]

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(glu_vals, glt_vals, biomass_vals)
ax.set_xlabel('glucose')
ax.set_ylabel('oxygen')
ax.set_zlabel('biomass')
ax.set_zlim3d(0,0.005)
plt.show()


# split into input (X) and output (y) variables
X = training_data[:,0:2]
y = training_data[:,2:5]
# define the keras model
model = Sequential()
model.add(Dense(4, input_shape=(2,), activation='relu'))
model.add(Dense(12, activation='relu'))
model.add(Dense(3, activation='relu'))
# compile the keras model
model.compile(loss='mse', optimizer='adam', metrics=['mae'])
# fit the keras model on the dataset
history = model.fit(X, y, epochs=250, batch_size=10)

plt.plot(history.history['mae'],'o',color='black')
plt.title('mean absolutue error')
plt.ylabel('mea')
plt.xlabel('number of epochs')
plt.show()



#%%
x_test_data = test_data[:,0:2]
y_test_data = test_data[:,2:5]


# evaluate the keras model
print('EVALUATION')
_, accuracy = model.evaluate(x_test_data, y_test_data)


#%%

sol = model.predict([[1000., 20]])
print(sol)

#save model
from keras2cpp import export_model
export_model(model, 'example.model')
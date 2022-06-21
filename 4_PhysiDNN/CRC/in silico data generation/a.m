target = 4988;
start = 2494;
now = start;

counter = 0;
while target > now
    now = now * (1 + (0.035)/60);
    counter = counter + 1;
end
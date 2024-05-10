#ifndef __M_PARTICLE_H__
#define __M_PARTICLE_H__

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

namespace MP
{
    float PARTICLE_SIZE = 1.0f;

    class Particle
    {

    public:
        // physics properties
        glm::vec3 position;
        glm::vec3 velocity;
        glm::vec3 acceleration;
        glm::vec3 force;
        float mass;
        float damping;
        // control properties
        float life;
        float age;
        bool alive;

    public:
        Particle() = default;
        // pos: position, vel: velocity, acc: acceleration, ma: mass, d: damping
        Particle(glm::vec3 pos, glm::vec3 vel, glm::vec3 acc, float ma, float d)
            : position(pos), velocity(vel), acceleration(acc), force(glm::vec3(0.0f)),
            mass(ma), damping(d), life(1.0f), age(0.0f), alive(true)
        {    }


        ~Particle()= default;

        void Update(){
            // todo: update physics
        };
    };


    void CreateParticle_Cube();
    // void CreateParticle_Sphere();
}

#endif // __M_PARTICLE_H__
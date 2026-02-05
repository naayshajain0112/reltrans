! Import FTorch library for interfacing with PyTorch
use ftorch, only : torch_model, torch_tensor, torch_kCUDA, torch_delete, &
                   torch_tensor_from_array, torch_model_load, &
                   torch_model_forward, torch_kCPU

implicit none

subroutine init_emulator()

    ! Generate an object to hold the Torch model
    type(torch_model) :: emulator
    
    ! Initialise the Torch model to be used
    torch_model_load(emulator, pathname_emulator, torch_kCUDA)

end subroutine init_emulator

subroutine run_model_inference(emulator, photarx, earx, nex, h, a, inc, rin, rout, gamma, logxi, &
                                Afe, lognep, Cutoff_s, boost)
    ! Set up array of n_inputs input tensors and array of n_outputs output tensors
    integer, parameter :: n_inputs = 10
    type(torch_tensor), dimension(n_inputs)  :: model_pars
    type(torch_tensor), dimension(nex)  :: model_egrid
    type(torch_tensor), dimension(nex) :: model_output_arr

    ! Set up the model inputs and output as Fortran arrays (pars, energy to evaluate, output)
    real, dimension(n_inputs), target    :: pars
    real, dimension(nex), target          :: egrid
    real, dimension(nex), target          :: output

    ! pass relevant parameters to fortran array before passing into emulator
    pars(1)  = h(1)         !height
    pars(2)  = a            !spin
    pars(3)  = inc          !inclination
    pars(4)  = rin          !inner disk radius
    pars(5)  = rout         !outer disk radius
    pars(6)  = gamma        !Gamma
    pars(7)  = logxi        !logxi
    pars(8)  = Afe          !Afe
    pars(9)  = lognep       !lognep
    pars(10) = Cutoff_s     !kTe

    ! Wrap Fortran data as no-copy Torch Tensors
    call torch_tensor_from_array(model_pars(1), pars, torch_kCUDA)
    call torch_tensor_from_array(model_egrid(1), earx, torch_kCUDA)
    call torch_tensor_from_array(model_output_arr(1), output, torch_kCPU)

    ! Run model forward method and infer reflection spectrum
    call torch_model_forward(emulator, model_pars, model_egrid, model_output_arr)

    ! Clean up
    call torch_delete(model_pars)
    call torch_delete(model_egrid)
    call torch_delete(model_output_arr)

    do i = 1, nex
        photarx(i) = output(i) * boost
    end do

end subroutine run_model_inference
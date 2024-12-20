//
// Created by Syl on 22/11/2024.
//

#ifndef FFT2_FFT2_IPP_H
#define FFT2_FFT2_IPP_H
#include <ipp/ippcore.h>
#include <ipp/ipps.h>
#include <memory>
#include <iostream>
#include <cassert>
namespace signalsmith { namespace fft2 {

namespace {

void checkStatus(IppStatus status) {
  assert(status == ippStsNoErr);
}

template<typename T>
struct State{};

template<>
struct State<float> {
public:
    explicit State(size_t size) {
        IppStatus status;
        const auto order = static_cast<int>(std::log2(size));
        Ipp8u* initBuffer{ nullptr };
        int fftSpecSize, fftInitBuffSize, fftWorkBuffSize;
        status = ippsFFTGetSize_C_32fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast, &fftSpecSize, &fftInitBuffSize, &fftWorkBuffSize);
        checkStatus(status);
        specBuffer = ippsMalloc_8u(fftSpecSize);
        spec = (IppsFFTSpec_C_32fc*)specBuffer;
        if(fftInitBuffSize != 0) {
            initBuffer = ippsMalloc_8u(fftInitBuffSize);
        }
        if(fftWorkBuffSize != 0) {
            workBuffer = ippsMalloc_8u(fftWorkBuffSize);
        }
        status = ippsFFTInit_C_32fc(&spec, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast, specBuffer, initBuffer);
        checkStatus(status);
        if(initBuffer) {
            ippsFree(initBuffer);
        }

    }

    ~State() noexcept {
        if(specBuffer) {
            ippsFree(specBuffer);
        }
        if(workBuffer) {
            ippsFree(workBuffer);
        }
    }

    IppsFFTSpec_C_32fc *spec{nullptr};
    Ipp8u *specBuffer{nullptr};
    Ipp8u *workBuffer{nullptr};

};

template<>
struct State<double> {
    explicit State(size_t size) {
        IppStatus status;
        const auto order = static_cast<int>(std::log2(size));
        Ipp8u* initBuffer{ nullptr };
        int fftSpecSize, fftInitBuffSize, fftWorkBuffSize;
        status = ippsFFTGetSize_C_64fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast, &fftSpecSize, &fftInitBuffSize, &fftWorkBuffSize);
        checkStatus(status);
        specBuffer = ippsMalloc_8u(fftSpecSize);
        spec = (IppsFFTSpec_C_64fc*)specBuffer;
        if(fftInitBuffSize != 0) {
            initBuffer = ippsMalloc_8u(fftInitBuffSize);
        }
        if(fftWorkBuffSize != 0) {
            workBuffer = ippsMalloc_8u(fftWorkBuffSize);
        }
        status = ippsFFTInit_C_64fc(&spec, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast, specBuffer, initBuffer);
        checkStatus(status);
        if(initBuffer) {
            ippsFree(initBuffer);
        }
    }

    ~State() noexcept {
        if(specBuffer) {
            ippsFree(specBuffer);
        }
        if(workBuffer) {
            ippsFree(workBuffer);
        }
    }

    IppsFFTSpec_C_64fc *spec{nullptr};
    Ipp8u *specBuffer{nullptr};
    Ipp8u *workBuffer{nullptr};

};
}

template<>
struct Pow2FFT<float> {
public:
    Pow2FFT(size_t size = 0) {
        if(size > 0) {
            resize(size);
        }
    }

    void resize(size_t size) {
        m_state = std::make_unique<State<float>>(size);
    }

    void fft(const std::complex<float>* input, std::complex<float>* output) {
        auto* ippComplexIn = reinterpret_cast<const Ipp32fc*>(input);
        auto* ippComplexOut = reinterpret_cast<Ipp32fc*>(output);
        const auto status = ippsFFTFwd_CToC_32fc(ippComplexIn, ippComplexOut, m_state->spec, m_state->workBuffer);
        checkStatus(status);
    }

    void ifft(const std::complex<float>* input, std::complex<float>* output) {
        auto* ippComplexIn = reinterpret_cast<const Ipp32fc*>(input);
        auto* ippComplexOut = reinterpret_cast<Ipp32fc*>(output);
        const auto status = ippsFFTInv_CToC_32fc(ippComplexIn, ippComplexOut, m_state->spec, m_state->workBuffer);
        checkStatus(status);
    }

private:
  std::unique_ptr<State<float>> m_state{ nullptr };
};

template<>
struct Pow2FFT<double> {
public:
    Pow2FFT(size_t size = 0) {
        if(size > 0) {
            resize(size);
        }
    }

    void resize(size_t size) {
        m_state = std::make_unique<State<double>>(size);
    }

    void fft(const std::complex<double>* input, std::complex<double>* output) {
        auto* ippComplexIn = reinterpret_cast<const Ipp64fc*>(input);
        auto* ippComplexOut = reinterpret_cast<Ipp64fc*>(output);
        const auto status = ippsFFTFwd_CToC_64fc(ippComplexIn, ippComplexOut, m_state->spec, m_state->workBuffer);
        checkStatus(status);
    }

    void ifft(const std::complex<double>* input, std::complex<double>* output) {
        auto* ippComplexIn = reinterpret_cast<const Ipp64fc*>(input);
        auto* ippComplexOut = reinterpret_cast<Ipp64fc*>(output);
        const auto status = ippsFFTInv_CToC_64fc(ippComplexIn, ippComplexOut, m_state->spec, m_state->workBuffer);
        checkStatus(status);
    }

private:
  std::unique_ptr<State<double>> m_state{ nullptr };
};

}}

#endif // FFT2_FFT2_IPP_H

import { onScopeDispose, ref } from 'vue'

/**
 * 复制反馈状态：showSuccess 置位 copied，duration 毫秒后自动复原；
 * showFailure 立即清除成功态（失败不显示成功，由调用方决定 UI 语义）。
 * timer 随作用域销毁自动清理，避免组件卸载后泄漏。
 */
export function useCopyFeedback(duration = 1200) {
  const copied = ref(false)
  let timer: ReturnType<typeof setTimeout> | null = null

  const clearTimer = () => {
    if (timer !== null) {
      clearTimeout(timer)
      timer = null
    }
  }

  const showSuccess = () => {
    clearTimer()
    copied.value = true
    timer = setTimeout(() => {
      copied.value = false
      timer = null
    }, duration)
  }

  const showFailure = () => {
    clearTimer()
    copied.value = false
  }

  onScopeDispose(clearTimer)

  return { copied, showSuccess, showFailure }
}

<template>
  <FloatingButton
    source-id="btt"
    :default-top="defaultTop"
    mode="match"
    :on-release="backToTop"
    :ariaLabel="t('backToTop')"
    :show="showBackToTop"
  >
    <i class="fas fa-arrow-up"></i>
  </FloatingButton>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import FloatingButton from '@/components/widgets/FloatingButton.vue'
import { scrollToTop } from '@/utils/scroll'

const { t } = useI18n()

const showBackToTop = ref(false)
const scrollObserver = ref<IntersectionObserver | null>(null)
const scrollSentinel = ref<HTMLDivElement | null>(null)

const defaultTop = (typeof window !== 'undefined' ? window.innerHeight : 1024) - 100

function backToTop() {
  scrollToTop()
}

function setupScrollObserver() {
  if (scrollObserver.value) return
  const sentinel = document.createElement('div')
  sentinel.style.cssText =
    'position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden'
  document.body.appendChild(sentinel)
  scrollSentinel.value = sentinel

  scrollObserver.value = new IntersectionObserver(
    ([entry]) => {
      showBackToTop.value = !entry.isIntersecting
    },
    { threshold: 0 },
  )
  scrollObserver.value.observe(sentinel)
}

onMounted(() => {
  setupScrollObserver()
})

onBeforeUnmount(() => {
  if (scrollObserver.value) {
    scrollObserver.value.disconnect()
    scrollObserver.value = null
  }
  if (scrollSentinel.value) {
    scrollSentinel.value.remove()
    scrollSentinel.value = null
  }
})
</script>

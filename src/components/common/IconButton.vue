<!-- 操作图标按钮：统一 stroke SVG 骨架 + 点击失焦（AGENTS.md 图标规范的单一实现） -->
<script setup lang="ts">
import { ICON_MARKS, type IconName } from '@/utils/icons'

withDefaults(defineProps<{ icon: IconName; ariaLabel: string; size?: number }>(), { size: 18 })

const blurFocus = (event: FocusEvent) => {
  ;(event.currentTarget as HTMLButtonElement | null)?.blur()
}
</script>

<template>
  <button type="button" class="icon-btn" :aria-label="ariaLabel" @focus="blurFocus">
    <svg
      class="icon-btn__svg"
      :width="size"
      :height="size"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      stroke-width="1.75"
      stroke-linecap="round"
      stroke-linejoin="round"
      aria-hidden="true"
      v-html="ICON_MARKS[icon]"
    />
  </button>
</template>

<style scoped>
.icon-btn__svg {
  display: block;
}
</style>

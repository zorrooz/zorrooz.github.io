<template>
  <div class="empty-state" v-reveal>
    <i v-if="icon" :class="['empty-state__icon', icon]" aria-hidden="true"></i>
    <p class="empty-state__text"><slot /></p>
    <router-link v-if="showHome" :to="homePath" class="empty-state__home">
      <i class="fas fa-house"></i>{{ t('backHome') }}
    </router-link>
  </div>
</template>

<script setup lang="ts">
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { toLocalePath } from '@/utils/navigation'

withDefaults(defineProps<{ icon?: string; showHome?: boolean }>(), {
  icon: '',
  showHome: false,
})

const { t } = useI18n()
const homePath = computed(() => toLocalePath('/'))
</script>

<style scoped>
.empty-state {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  gap: var(--sp-3);
  padding: var(--sp-16) var(--sp-6);
  text-align: center;
  background: var(--surface-2);
  border: 1px solid var(--line);
  border-radius: var(--radius);
}

.empty-state__icon {
  font-size: 26px;
  color: var(--fg-3);
}

.empty-state__text {
  margin: 0;
  font-size: var(--text-md);
  color: var(--fg-2);
}

.empty-state__home {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  margin-top: var(--sp-2);
  padding: 8px 18px;
  font-size: var(--text-sm);
  font-weight: 500;
  color: var(--primary);
  background: transparent;
  border: 1px solid var(--line);
  border-radius: var(--radius-btn);
  text-decoration: none;
  transition:
    color 0.14s ease,
    border-color 0.14s ease,
    background-color 0.14s ease;
}

.empty-state__home:hover {
  color: var(--primary);
  border-color: color-mix(in srgb, var(--primary) 35%, transparent);
  background: var(--tint);
}
</style>
